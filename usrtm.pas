unit usrtm;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, ContNrs, Zipper, SyncObjs;

type
  TSRTMSegment = class
  private
    FStream: TMemoryStream;
    FBlockSize: integer;
    FLat: integer;
    FLon: integer;

    procedure DoCreateOutZipStream(Sender: TObject; var AStream: TStream;
    {%H-}AItem: TFullZipFileEntry);
    procedure DoDoneOutZipStream(Sender: TObject; var AStream: TStream;
    {%H-}AItem: TFullZipFileEntry);

  public
    constructor Create(const ALat, ALon: double; const AHDMode: boolean;
      const ACacheDir: string = '');
    destructor Destroy; override;

    function GetElevation(const ALat, ALon: double): integer;
    function HasElevation(const ALat, ALon: double): boolean;
  end;

  TSRTM = class(TThread)
  private
    FSegments: TObjectList;
    FCacheDir: string;
    FLat, FLon: double;
    FDestLat, FDestLon: double;

    FDistance: double;
    FIsLoading: TEvent;
    FCanRead: TEvent;
    FUseRadius: boolean;

  protected
    procedure Execute; override;
  public
    constructor Create(const ACacheDir: string = '');
    destructor Destroy; override;

    function Load(const ALat, ALon: double; const ADistanceMiles: double = 50): boolean;
      overload;
    function Load(const ASourceLat, ASourceLon, ADestLat, ADestLon: double): boolean;
      overload;

    function GetElevation(const ALat, ALon: double): integer;
  end;

implementation

uses fphttpclient, LazFileUtils;

const
  EARTHRADIUS_MILES = 3959.0;
  DEG2RAD = 1.74532925199e-02;

  SRTM1_URL = 'http://srtm.kurviger.de/SRTM1/';
  SRTM3_URL = 'http://srtm.kurviger.de/SRTM3/';

function GetHGTFile(ALat, ALon: double): string;
var
  sLat, sLon: char;
  lonOffset: byte;
begin
  if Alat >= 0 then
    sLat := 'N'
  else
    sLat := 'S';

  if Alon >= 0 then
  begin
    sLon := 'E';
    lonOffset := 0;
  end
  else
  begin
    sLon := 'W';
    lonOffset := 1;
  end;

  Result := format('%s%.02d%s%.03d.hgt', [sLat, trunc(abs(Alat)),
    sLon, trunc(abs(Alon) + lonOffset)]);
end;

function GetUSRegion(const ABaseUrl: string; const ALat, ALon: double): string;
  //Return the SRTM1 region number of a given lat, lon.
  //http://dds.cr.usgs.gov/srtm/version2_1/SRTM1/Region_definition.jpg
var
  region: integer;
begin
  if (38 <= Alat) and (Alat < 50) and (-125 <= Alon) and (Alon < -111) then
    region := 1
  else if (38 <= Alat) and (Alat < 50) and (-111 <= Alon) and (Alon < -97) then
    region := 2
  else if (38 <= Alat) and (Alat < 50) and (-97 <= Alon) and (Alon < -83) then
    region := 3
  else if (28 <= Alat) and (Alat < 38) and (-123 <= Alon) and (Alon < -100) then
    region := 4
  else if (25 <= Alat) and (Alat < 38) and (-100 <= Alon) and (Alon < -83) then
    region := 5
  else if (17 <= Alat) and (Alat < 48) and (-83 <= Alon) and (Alon < -64) then
    region := 6
  else if (-15 <= Alat) and (Alat < 60) and (((172 <= Alon) and (Alon < 180)) or
    ((-180 <= Alon) and (Alon < -129))) then
    region := 7
  else
    region := 0; //unknown region

  if region <> 0 then
    Result := ABaseUrl + format('Region_%0.2d/', [region])
  else
    Result := '';
end;

function GetWorldRegion(const ABaseUrl: string; const ALat, ALon: double): string;
const
  Regions: array[0..5] of string =
    ('Africa', 'Australia', 'Eurasia', 'Islands', 'North_America', 'South_America');
{$i srtm_regions.inc}
var
  str, helper: string;
  i: integer;
begin
  str := StringReplace(GetHGTFile(ALat, ALon), '.hgt', ',', []);
  for i := 0 to 5 do
  begin
    case i of
      0: helper := cAfrica;
      1: helper := cAustralia;
      2: helper := cEurasia;
      3: helper := cIslands;
      4: helper := cNorthAmerica;
      else
        helper := cSouthAmerica;
    end;
    if pos(str, helper) > 0 then
    begin
      Result := ABaseUrl + Regions[i] + '/';
      exit;
    end;
  end;
  Result := '';
end;

function IsValidFile(const AFilename: string): boolean;
var
  Stream: TFileStream;
begin
  Result := FileExists(AFilename);
  if Result then
  begin
    Stream := TFileStream.Create(AFilename, fmOpenRead);
    Result := Stream.Size > 0;
    Stream.Free;
  end;
end;

{ TSRTMSegment }
constructor TSRTMSegment.Create(const ALat, ALon: double; const AHDMode: boolean;
  const ACacheDir: string = '');
var
  Stream: TMemoryStream;
  http: TFPHttpClient;
  Filename: string;
  DownloadUrl: string;
  ZipFile: TUnZipper;
  CacheDir: string;
begin
  CacheDir := ACacheDir;
  if CacheDir = '' then
    CacheDir := IncludeTrailingPathDelimiter(GetTempDir);

  try
    if not DirectoryExists(CacheDir) then
      ForceDirectory(CacheDir);
  except
  end;

  FLat := trunc(ALat);
  FLon := trunc(ALon);

  Filename := GetHgtFile(ALat, ALon) + '.zip';
  FStream := TMemoryStream.Create;

  if AHDMode then
    DownloadUrl := GetUSRegion(SRTM1_URL, ALat, ALon)
  else
    DownloadUrl := '';

  if DownloadUrl = '' then  //Get SD Mode
    DownloadUrl := GetWorldRegion(SRTM3_URL, ALat, ALon);


  if not FileExists(CacheDir + Filename) then
  begin
    Stream := TMemoryStream.Create;
    http := TFPHttpClient.Create(nil);
    try
      try
        http.Get(DownloadUrl + Filename, Stream);
      except
      end;
      if http.ResponseStatusCode = 404 then //typo in the sources, where a . is missing
        try
          http.Get(DownloadUrl + StringReplace(Filename, '.hgt.zip',
            'hgt.zip', []), Stream);
        except
        end;
      if Stream.Size > 0 then
        Stream.SaveToFile(CacheDir + Filename);
    finally
      Stream.Free;
      http.Free;
    end;
  end;

  if IsValidFile(CacheDir + Filename) then
  begin
    ZipFile := TUnZipper.Create;
    try
      ZipFile.FileName := CacheDir + Filename;
      ZipFile.OnCreateStream := @DoCreateOutZipStream;
      ZipFile.OnDoneStream := @DoDoneOutZipStream;
      ZipFile.UnZipAllFiles;
    finally
      ZipFile.Free;
    end;
  end;

end;

destructor TSRTMSegment.Destroy;
begin
  FreeAndNil(FStream);
  inherited;
end;

procedure TSRTMSegment.DoCreateOutZipStream(Sender: TObject;
  var AStream: TStream; AItem: TFullZipFileEntry);
begin
  FBlockSize := 0;
  AStream := FStream;
end;

procedure TSRTMSegment.DoDoneOutZipStream(Sender: TObject; var AStream: TStream;
  AItem: TFullZipFileEntry);
begin
  AStream.Position := 0;
  if AStream.Size = 1201 * 1201 * 2 then // SRTM-3
    FBlockSize := 1201
  else
  if AStream.Size = 3601 * 3601 * 2 then // SRTM-1
    FBlockSize := 3601
  else
    FBlockSize := 0;
end;

function TSRTMSegment.HasElevation(const ALat, ALon: double): boolean;
begin
  Result := (Trunc(ALat) = FLat) and (Trunc(ALon) = FLon);
end;

function TSRTMSegment.GetElevation(const ALat, ALon: double): integer;
var
  lat_row, lon_row: integer;
  Data: array[0..1] of byte;
begin
  Result := 0;
  if (not HasElevation(ALat, ALon)) then
    exit;

  if (not assigned(FStream)) or (FBlockSize = 0) then
    exit;

  Data[0] := 0;
  Data[1] := 0;

  // Calculate location in the file.
  lat_row := round((ALat - Trunc(ALat)) * (FBlockSize - 1));
  lon_row := round((ALon - Trunc(ALon)) * (FBlockSize - 1));

  if lat_row < 0 then
    lat_row := abs(lat_row)
  else
    lat_row := FBlockSize - 1 - lat_row;

  FStream.Position := (int64(lat_row) * FBlockSize + lon_row) * 2;
  FStream.Read(Data[0], 2);

  Result := Data[0] shl 8 or Data[1];

  if (Data[0] and 128 > 0) then
    Result := Result - $10000;
  if (Result < -32768) then
    Result := Result - 32768;
  if (Result > 32767) then
    Result := 32767;
end;


{ TSRTM }
constructor TSRTM.Create(const ACacheDir: string = '');
begin
  inherited Create(False);
  FSegments := TObjectList.Create(True);
  FCacheDir := ACacheDir;
  FIsLoading := TEvent.Create(nil, True, False, '');
  FCanRead := TEvent.Create(nil, True, False, '');
end;

destructor TSRTM.Destroy;
begin
  Terminate;
  FIsLoading.SetEvent;
  WaitFor;

  FIsLoading.Free;
  FCanRead.Free;
  FSegments.Free;
  inherited;
end;

function TSRTM.Load(const ASourceLat, ASourceLon, ADestLat, ADestLon: double): boolean;
begin
  if FIsLoading.WaitFor(0) = TWaitResult.wrSignaled then
  begin
    Result := False;
    exit;
  end;


  if (not FUseRadius) and (FLat = ASourceLat) and (FLon = ASourceLon) and
    (FDestLat = ADestLat) and (ADestLon = ADestLon) then
  begin
    // do nothing
  end
  else
  begin
    FLat := ASourceLat;
    FLon := ASourceLon;
    FDestLat := ADestLat;
    FDestLon := ADestLon;
    FUseRadius := False;
    FDistance := 0;
    FIsLoading.SetEvent;
  end;

  Result := True;
end;

function TSRTM.Load(const ALat, ALon: double;
  const ADistanceMiles: double = 50): boolean;
begin
  if FIsLoading.WaitFor(0) = TWaitResult.wrSignaled then
  begin
    Result := False;
    exit;
  end;

  if (FUseRadius) and (FLat = ALat) and (FLon = ALon) and
    (FDistance = ADistanceMiles) then
  begin
    // do nothing
  end
  else
  begin
    FLat := ALat;
    FLon := ALon;
    FUseRadius := True;
    FDistance := ADistanceMiles;
    FIsLoading.SetEvent;
  end;
  Result := True;
end;

procedure TSRTM.Execute;
var
  east, north, south, west: double;
  ilat, ilon: integer;
begin
  while not terminated do
  begin
    FIsLoading.WaitFor(Infinite);
    FCanRead.ResetEvent;

    if Terminated then
      exit;

    FSegments.Clear;
    if FUseRadius then
    begin
      north := FLat - (FDistance / EARTHRADIUS_MILES) / DEG2RAD;
      south := FLat + (FDistance / EARTHRADIUS_MILES) / DEG2RAD;
      east := FLon - (FDistance / EARTHRADIUS_MILES) / DEG2RAD / Cos(FLat * DEG2RAD);
      west := FLon + (FDistance / EARTHRADIUS_MILES) / DEG2RAD / Cos(FLat * DEG2RAD);
    end
    else
    begin
      if FLat < FDestLat then
        north := FLat
      else
        north := FDestLat;
      if FLat > FDestLat then
        south := FLat
      else
        south := FDestLat;
      if FLon < FDestLon then
        east := FLon
      else
        east := FDestLon;
      if FLon > FDestLon then
        west := FLon
      else
        west := FDestLon;
    end;
    for ilat := trunc(north) to trunc(south) do
      for ilon := trunc(east) to trunc(west) do
        FSegments.Add(TSRTMSegment.Create(ilat, ilon, True, FCacheDir));

    FCanRead.SetEvent;
    FIsLoading.ResetEvent;
  end;
end;

function TSRTM.GetElevation(const ALat, ALon: double): integer;
var
  i: integer;
begin
  FCanRead.WaitFor(Infinite);

  for i := 0 to FSegments.Count - 1 do
    if TSRTMSegment(FSegments[i]).HasElevation(ALat, ALon) then
    begin
      Result := TSRTMSegment(FSegments[i]).GetElevation(ALat, ALon);
      exit;
    end;
  Result := 0;
end;

end.
