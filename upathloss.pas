unit upathloss;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, ContNrs, SyncObjs, uitwom, usrtm;

type
  TOnGetElevation = procedure(Sender: TObject; ALat, ALon: double;
    out Elevation: double) of object;

  TSite = record
    Lat, Lon: double;
    Alt: double;
    Caption: string;
  end;

  TPathItem = class
  private
    FLat: double;
    FLon: double;
    FElevation: double;
    FDistance: double;
  public
    constructor Create(const ALat, ALon, AElevation, ADistance: double);
    function AsString: string;

    property Lat: double read FLat;
    property Lon: double read FLon;
    property Elevation: double read FElevation write FElevation;
    property Distance: double read FDistance;
  end;

  TPath = class(TObjectList)
  public
    constructor Create(const ASource, ADestination: TSite;
      const APixelPerDegree: integer; const AGetElevation: TOnGetElevation);
  end;

  TAntennaPattern = array[0..360, 0..1000] of double;

  TSettings = class
  private
    Feps_dielect: double;
    Fsgm_conductivity: double;
    Feno_ns_surfref: double;
    Ffrq_mhz: double;
    Fconf: double;
    Frel: double;
    Ferp: double;
    Fradio_climate: integer;
    Fpol: integer;
    Fantenna_pattern: TAntennaPattern;

    FEnvironment: TEnvironment; //HATA
  public
    constructor Create;
    procedure LoadAntennaPattern(const AAzimute, AElevation: TStrings);

    property eps_dielect: double read Feps_dielect write Feps_dielect;
    property sgm_conductivity: double read Fsgm_conductivity write Fsgm_conductivity;
    property eno_ns_surfref: double read Feno_ns_surfref write Feno_ns_surfref;
    property frq_mhz: double read Ffrq_mhz write Ffrq_mhz;
    property conf: double read Fconf write Fconf;
    property rel: double read Frel write Frel;
    property erp: double read Ferp write Ferp;
    property radio_climate: integer read Fradio_climate write Fradio_climate;
    property pol: integer read Fpol write FPol;

    property Environment: TEnvironment read FEnvironment write FEnvironment;

    property antenna_pattern: TAntennaPattern read Fantenna_pattern;
  end;

  TPathlossModel = (pmLongleyRice, pmHATA, pmECC33, pmSUI,
    pmCost231, pmFreespace, pmITWOM3, pmEricsson, pmPlainEarth, pmEgli, pmSoil);

  TRangeSettings = record
    min_west, max_west: double;
    min_north, max_north: double;
    eastwest: boolean;
    Source: TSite;
    Altitude: double; // Antenna height
    los: boolean;
    mask_value: byte;
  end;

  T2DArray = array of byte;

  TPathloss = class
  private
    FLock: TCriticalSection;
    FComplete: TEvent;
    FAbort: TEvent;

    FTaskCount: integer;

    FOnComplete: TNotifyEvent;

    FMaskValue: byte;
    FClutter: double;
    FMaxRange: double;

    FMask: T2DArray;
    FSignal: T2DArray;
    FAltitude: array of integer;

    FModel: TPathlossModel;
    FUseDBm: boolean;
    FGotElevationPattern: boolean;

    FSettings: TSettings;

    FHottest: byte;

    FMinWest: double;
    FMaxWest: double;
    FMinNorth: double;
    FMaxNorth: double;

    FHeight: integer;
    FHDMode: boolean;

    FMinElevation: double;
    FMaxElevation: double;

    FSource: TSite;

    procedure SetHDMode(const AValue: boolean);
    procedure SetMaxRange(const AValue: double);
    procedure Propagate(const ARange: TRangeSettings);

    function PutMask(const ALat, ALon: double; const AValue: integer): integer;
    function OrMask(const ALat, ALon: double; const AValue: integer): integer;
    function GetMask(const ALat, ALon: double): integer;

    procedure PutSignal(const ALat, ALon: double; const ASignal: byte); virtual;
    function GetSignal(const ALat, ALon: double): byte;

    procedure PlotLOSPath(const ASource, ADestination: TSite; const AMaskValue: byte);
    procedure PlotPropPath(const ASource, ADestination: TSite; const AMaskValue: byte);

    function GetPixelPerDegree: integer;
    procedure CompleteTask();
  public
    constructor Create;
    destructor Destroy; override;

    procedure Terminate;
    procedure WaitFor(const ATimeOut: cardinal = INFINITE);

    function IsCalulating: boolean;

    function SiteToOffset(const ALat, ALon: double;
      out X, Y, Offset: integer): boolean; overload;
    function SiteToOffset(const ALat, ALon: double; out Offset: integer): boolean;
      overload;

    function GetMaskValue(const AOffset: integer): byte;
    function GetSignalValue(const AOffset: integer): byte;

    procedure GetElevation(Sender: TObject; ALat, ALon: double;
      out Elevation: double);
    procedure ReadAltitudeMap(const ASource: TSite);
    function GetAverageHeight(const AX, AY: integer; ASize: integer = 2): integer;

    procedure Calculate(const ASource: TSite; const AAltitude: double;
      const ALineOfSight: boolean = False);

    property Clutter: double read FClutter write FClutter;
    property GotElevationPattern: boolean read FGotElevationPattern
      write FGotElevationPattern;

    property Height: integer read FHeight;
    property HDMode: boolean read FHDMode write SetHDMode;
    property Model: TPathlossModel read FModel write FModel;

    property PixelPerDegree: integer read GetPixelPerDegree;
    property Settings: TSettings read FSettings;
    property UseDBm: boolean read FUseDBm write FUseDBM;

    property MaxRange: double read FMaxRange write SetMaxRange;
    property MinNorth: double read FMinNorth;
    property MaxNorth: double read FMaxNorth;
    property MinWest: double read FMinWest;
    property MaxWest: double read FMaxWest;
    property MaxElevation: double read FMaxElevation;
    property MinElevation: double read FMinElevation;

    property OnComplete: TNotifyEvent read FOnComplete write FOnComplete;
  end;

  TCalculatorItem = class
  private
    FLat: double;
    FLon: double;
    FSignal: byte;
  public
    constructor Create(const ALat, ALon: double; const ASignal: byte);
    property Lat: double read FLat;
    property Lon: double read FLon;
    property Signal: byte read FSignal;
  end;

  TCalculator = class(TObjectList)
  protected
    FHDMode: boolean;
    FClutter: double;
    FSRTM: TSRTM;
    procedure GetElevation(Sender: TObject; ALat, ALon: double;
      out Elevation: double);
    function GetPixelPerDegree: integer;
  public
    constructor Create;
    destructor Destroy; override;

    procedure Calculate(const ASource, ADestination: TSite); virtual; abstract;

    property HDMode: boolean read FHDMode write FHDMode;
    property Clutter: double read FClutter write FClutter;
  end;

  TLOSCalculator = class(TCalculator)
  public
    procedure Calculate(const ASource, ADestination: TSite); override;
  end;

  TPropCalculator = class(TCalculator)
  private
    FGotElevationPattern: boolean;
    FModel: TPathlossModel;
    FUseDBm: boolean;
    FSettings: TSettings;
  public
    constructor Create;
    destructor Destroy; override;

    procedure Calculate(const ASource, ADestination: TSite); override;

    property GotElevationPattern: boolean read FGotElevationPattern
      write FGotElevationPattern;

    property Model: TPathLossModel read FModel write FModel;
    property UseDBm: boolean read FUseDBm write FUseDBm;
    property Settings: TSettings read FSettings;
  end;

implementation

uses Math, LazFileUtils;

const
  DEG2RAD = 1.74532925199e-02;
  TWOPI = 6.283185307179586;
  HALFPI = 1.570796326794896;
  FOUR_THIRDS = 1.3333333333333;
  FEET_PER_MILE = 5280.0;
  METERS_PER_FOOT = 0.3048;
  METERS_PER_MILE = 1609.344;

  EARTHRADIUS = 20902230.97;
  EARTHRADIUS_MILES = 3959.0;
  RENDER_AROUND_RANGE = 5;

type
  TProcessThread = class(TThread)
  protected
    FParent: TPathLoss;
    FRange: TRangeSettings;

    procedure Execute; override;
  public
    constructor Create(const ARange: TRangeSettings; const AParent: TPathloss);
    destructor Destroy; override;
  end;

constructor TProcessThread.Create(const ARange: TRangeSettings;
  const AParent: TPathloss);
begin
  inherited Create(False);
  FParent := APArent;
  FRange := Arange;
  FreeOnTerminate := True;
end;

procedure TProcessThread.Execute;
begin
  if assigned(FParent) then
    FParent.Propagate(FRange);
end;

destructor TPRocessThread.Destroy;
begin
  if assigned(FParent) then
    FParent.CompleteTask();
  inherited;
end;

function GetToken(var Line: string; ch: char = #9): string;
var
  i: integer;
begin
  i := pos(ch, Line);
  if i > 0 then
  begin
    Result := copy(Line, 1, i - 1);
    Delete(Line, 1, i);
  end
  else
  begin
    Result := Line;
    Line := '';
  end;
end;


function GetDistance(const ASite1, ASite2: TSite): double;
var
  lat1, lon1, lat2, lon2: double;
begin
    (* This function returns the great circle distance
       in miles between any two site locations. *)

  lat1 := ASite1.lat * DEG2RAD;
  lon1 := ASite1.lon * DEG2RAD;
  lat2 := ASite2.lat * DEG2RAD;
  lon2 := ASite2.lon * DEG2RAD;

  Result := EARTHRADIUS_MILES * arccos(sin(lat1) * sin(lat2) +
    cos(lat1) * cos(lat2) * cos((lon1) - (lon2)));
end;

function FindPointAtDistanceFrom(const AStartPoint: TSite;
  const ADistance, ABearing: double): TSite; //in km
const
  radius = 6378.14;
var
  lat1, lon1, brng, lat2, lon2: double;
begin
  lat1 := AStartPoint.Lat * DEG2RAD;
  lon1 := AStartPoint.Lon * DEG2RAD;
  brng := ABearing * DEG2RAD;
  lat2 := Arcsin(Sin(lat1) * Cos(ADistance / radius) + Cos(lat1) *
    Sin(ADistance / radius) * Cos(brng));
  lon2 := lon1 + Arctan2(Sin(brng) * Sin(ADistance / radius) *
    Cos(lat1), Cos(ADistance / radius) - Sin(lat1) * Sin(lat2));
  Result.lat := lat2 / DEG2RAD;
  Result.lon := lon2 / DEG2RAD;
end;

function GetAzimuth(const ASource, ADestination: TSite): double;
var
  dest_lat, dest_lon, src_lat, src_lon, beta, diff, num, den, fraction: double;
begin
    (* This function returns the azimuth (in degrees) to the
       destination as seen from the location of the source. *)
  dest_lat := ADestination.lat * DEG2RAD;
  dest_lon := ADestination.lon * DEG2RAD;

  src_lat := ASource.lat * DEG2RAD;
  src_lon := ASource.lon * DEG2RAD;

  (* Calculate Surface Distance *)

  beta := arccos(sin(src_lat) * sin(dest_lat) + cos(src_lat) *
    cos(dest_lat) * cos(src_lon - dest_lon));

  (* Calculate Azimuth *)

  num := sin(dest_lat) - (sin(src_lat) * cos(beta));
  den := cos(src_lat) * sin(beta);
  if den = 0 then
    den := 1;
  fraction := num / den;

  (* Trap potential problems in acos() due to rounding *)

  if (fraction >= 1.0) then
    fraction := 1.0
  else if (fraction <= -1.0) then
    fraction := -1.0;

  (* Calculate azimuth *)

  Result := arccos(fraction);

  (* Reference it to True North *)

  diff := dest_lon - src_lon;

  if (diff <= -PI) then
    diff := diff + TWOPI;
  if (diff >= PI) then
    diff := diff - TWOPI;

  if (diff > 0.0) then
    Result := TWOPI - Result;

  Result := (Result / DEG2RAD);
end;

function arccos2(const x, y: double): double;
begin
  (* This function implements the arc cosine function,
     returning a value between 0 and TWOPI. *)
  Result := 0;
  if (y > 0.0) then
    Result := arccos(x / y)
  else if (y < 0.0) then
    Result := PI + arccos(x / y);
end;

function LonDiff(const lon1, lon2: double): double;
begin
  (* This function returns the short path longitudinal
     difference between longitude1 and longitude2
     as an angle between -180.0 and +180.0 degrees.
     If lon1 is west of lon2, the result is positive.
     If lon1 is east of lon2, the result is negative. *)

  Result := lon1 - lon2;

  if (Result <= -180.0) then
    Result := Result + 360.0;

  if (Result >= 180.0) then
    Result := Result - 360.0;
end;

{ TPathItem }
constructor TPathItem.Create(const ALat, ALon, AElevation, ADistance: double);
begin
  FLat := ALat;
  FLon := ALon;
  FElevation := AElevation;
  FDistance := ADistance;
end;

function TPathItem.AsString: string;
begin
  Result := format('Latitude: %0.2f'#9'Longitude: %0.2f'#9 +
    'Elevation: %0.2f'#9'Distance: %0.2f', [FLat, FLon, FElevation, FDistance]);
end;

{ TPath }
constructor TPath.Create(const ASource, ADestination: TSite;
  const APixelPerDegree: integer; const AGetElevation: TOnGetElevation);
var
  azimuth, distance, lat1, lon1, beta, den, num, lat2, lon2, total_distance,
  dx, dy, path_length, miles_per_sample, samples_per_radian: double;
  elevation: double;
begin
  inherited Create(True);

  (* This function generates a sequence of latitude and
       longitude positions between source and destination
       locations along a great circle path, and stores
       elevation and distance information for points
       along that path in the "path" structure. *)

  samples_per_radian := 68755.0;

  lat1 := Asource.lat * DEG2RAD;
  lon1 := Asource.lon * DEG2RAD;

  lat2 := Adestination.lat * DEG2RAD;
  lon2 := Adestination.lon * DEG2RAD;
  Elevation := 0;
  samples_per_radian := APixelPerDegree * 57.295833;

  azimuth := GetAzimuth(ASource, ADestination) * DEG2RAD;

  total_distance := GetDistance(ASource, ADestination);

  if (total_distance > (30.0 / APixelPerDegree)) then (* > 0.5 pixel distance *)
  begin
    dx := samples_per_radian * arccos(cos(lon1 - lon2));
    dy := samples_per_radian * arccos(cos(lat1 - lat2));

    path_length := sqrt((dx * dx) + (dy * dy));    (* Total number of samples *)

    miles_per_sample := total_distance / path_length;  (* Miles per sample *)
  end
  else
  begin
    dx := 0.0;
    dy := 0.0;
    path_length := 0.0;
    miles_per_sample := 0.0;
    total_distance := 0.0;

    lat1 := lat1 / DEG2RAD;
    lon1 := lon1 / DEG2RAD;

    if assigned(AGetElevation) then
      AGetElevation(self, ASource.lat, ASource.Lon, Elevation);
    Add(TPathItem.Create(lat1, lon1, Elevation, 0.0));
  end;

  distance := 0.0;
  while (total_distance <> 0.0) and (distance <= total_distance) do
  begin
    distance := miles_per_sample * Count;
    beta := distance / EARTHRADIUS_MILES;
    lat2 := arcsin(sin(lat1) * cos(beta) + cos(azimuth) * sin(beta) * cos(lat1));
    num := cos(beta) - (sin(lat1) * sin(lat2));
    den := cos(lat1) * cos(lat2);


    if (azimuth = 0.0) and (beta > HALFPI - lat1) then
      lon2 := lon1 + PI
    else if (azimuth = HALFPI) and (beta > HALFPI + lat1) then
      lon2 := lon1 + PI
    else if (abs(num / den) > 1.0) then
      lon2 := lon1
    else
    begin
      if ((PI - azimuth) >= 0.0) then
        lon2 := lon1 - arccos2(num, den)
      else
        lon2 := lon1 + arccos2(num, den);
    end;

    lat2 := lat2 / DEG2RAD;
    lon2 := lon2 / DEG2RAD;

    if assigned(AGetElevation) then
      AGetElevation(self, lat2, lon2, Elevation);

    Add(TPathItem.Create(lat2, lon2, Elevation, distance));
  end;

  (* Make sure exact destination point is recorded at path.length-1 *)
  if assigned(AGetElevation) then
    AGetElevation(self, ADestination.lat, ADestination.lon, Elevation);
  Add(TPathItem.Create(Adestination.lat, Adestination.lon, Elevation,
    total_distance));
end;


{ TSettings }
constructor TSettings.Create;
begin
  Feps_dielect := 15.0;  // Farmland
  Fsgm_conductivity := 0.005;  // Farmland
  Feno_ns_surfref := 301.0;
  Ffrq_mhz := 200.0;  // Deliberately too low
  Fradio_climate := 5;  // continental
  Fpol := 1;    // vert
  Fconf := 0.50;
  Frel := 0.50;
  Ferp := 0.0;    // will default to Path Loss
  fillchar(Fantenna_pattern, sizeof(Fantenna_pattern), 0);
end;

procedure TSettings.LoadAntennaPattern(const AAzimute, AElevation: TStrings);
var
  w, i, x, y: integer;
  line: string;
  rotation: double;
  amplitude: double;
  elevation: double;
  slant_angle, azimuth_pattern, azimuth: array[0..360] of double;
  el_pattern: array[0..10000] of double;
  read_count: array[0..10000] of integer;
  elevation_pattern: array[0..360, 0..1000] of double;
  tilt_increment: double;
  xx, mechanical_tilt, tilt_azimuth: double;

  last_index, next_index: integer;
  valid1, valid2: double;
  span, delta: double;
  tilt: double;

  az, sum: double;
  a, b, z: integer;
  got_elevation_pattern, got_azimuth_pattern: boolean;
begin
  got_azimuth_pattern := False;
  got_elevation_pattern := False;

  if assigned(AAzimute) and (AAzimute.Count > 0) then
  begin
    rotation := StrToFloat(AAzimute[0]);
    for i := 0 to high(azimuth) do
    begin
      azimuth[i] := 0;
      read_count[i] := 0;
    end;
    for i := 1 to AAzimute.Count - 1 do
    begin
      line := AAzimute[i];
      x := round(StrToFloat(gettoken(line, #9)));
      amplitude := StrToFloat(gettoken(line, #9));
      if (x >= 0) and (x <= 360) then
      begin
        azimuth[x] := azimuth[x] + amplitude;
        Inc(read_count[x], 1);
      end;
    end;

    (* Handle 0=360 degree ambiguity *)

    if ((read_count[0] = 0) and (read_count[360] <> 0)) then
    begin
      read_count[0] := read_count[360];
      azimuth[0] := azimuth[360];
    end;

    if ((read_count[0] <> 0) and (read_count[360] = 0)) then
    begin
      read_count[360] := read_count[0];
      azimuth[360] := azimuth[0];
    end;

    (* Average pattern values in case more than
       one was read for each degree of azimuth. *)

    for x := 0 to 360 do
      if (read_count[x] > 1) then
        azimuth[x] := azimuth[x] / read_count[x];

    (* Interpolate missing azimuths
       to completely fill the array *)

    last_index := -1;
    next_index := -1;

    for x := 0 to 360 do
    begin
      if (read_count[x] <> 0) then
      begin
        if (last_index = -1) then
          last_index := x
        else
          next_index := x;
      end;

      if (last_index <> -1) and (next_index <> -1) then
      begin
        valid1 := azimuth[last_index];
        valid2 := azimuth[next_index];

        span := next_index - last_index;
        delta := (valid2 - valid1) / span;

        for y := last_index + 1 to next_index - 1 do
          azimuth[y] := azimuth[y - 1] + delta;

        last_index := y;
        next_index := -1;
      end;
    end;

    (* Perform azimuth pattern rotation
       and load azimuth_pattern[361] with
       azimuth pattern data in its final form. *)

    for x := 0 to 359 do
    begin
      y := x + round(rotation);

      if (y >= 360) then
        y := y - 360;

      azimuth_pattern[y] := azimuth[x];
    end;

    azimuth_pattern[360] := azimuth_pattern[0];

    got_azimuth_pattern := True;
  end;

  if assigned(AElevation) and (AElevation.Count > 0) then
  begin
    (* Clear azimuth pattern array *)

    for x := 0 to 10000 do
    begin
      el_pattern[x] := 0.0;
      read_count[x] := 0;
    end;
    line := AElevation[0];
    mechanical_tilt := round(StrToFloat(gettoken(line, #9)));
    tilt_azimuth := round(StrToFloat(gettoken(line, #9)));
    for i := 0 to AElevation.Count - 1 do
    begin
      line := AElevation[i];
      elevation := round(StrToFloat(gettoken(line, #9)));
      amplitude := round(StrToFloat(gettoken(line, #9)));

      x := round(100.0 * (elevation + 10.0));

      if (x >= 0) and (x <= 10000) then
      begin
        el_pattern[x] := el_pattern[x] + amplitude;
        Inc(read_count[x]);
      end;
    end;

    (* Average the field values in case more than
       one was read for each 0.01 degrees of elevation. *)

    for x := 0 to 10000 do
      if (read_count[x] > 1) then
        el_pattern[x] := el_pattern[x] / read_count[x];

    (* Interpolate between missing elevations (if
       any) to completely fill the array and provide
       radiated field values for every 0.01 degrees of
       elevation. *)

    last_index := -1;
    next_index := -1;

    for x := 0 to 10000 do
    begin
      if (read_count[x] <> 0) then
      begin
        if (last_index = -1) then
          last_index := x
        else
          next_index := x;
      end;

      if (last_index <> -1) and (next_index <> -1) then
      begin
        valid1 := el_pattern[last_index];
        valid2 := el_pattern[next_index];

        span := next_index - last_index;
        delta := (valid2 - valid1) / span;

        for y := last_index + 1 to next_index - 1 do
          el_pattern[y] :=
            el_pattern[y - 1] + delta;

        last_index := y;
        next_index := -1;
      end;
    end;
      (* Fill slant_angle[] array with offset angles based
       on the antenna's mechanical beam tilt (if any)
       and tilt direction (azimuth). *)

    if (mechanical_tilt = 0.0) then
    begin
      for x := 0 to 360 do
        slant_angle[x] := 0.0;
    end
    else
    begin
      tilt_increment := mechanical_tilt / 90.0;

      for x := 0 to 360 do
      begin
        xx := x;
        y := round(tilt_azimuth + xx);

        while (y >= 360) do
          Dec(y, 360);

        while (y < 0) do
          Inc(y, 360);

        if (x <= 180) then
          slant_angle[y] := -(tilt_increment * (90.0 - xx));

        if (x > 180) then
          slant_angle[y] := -(tilt_increment * (xx - 270.0));
      end;
    end;

    slant_angle[360] := slant_angle[0];  (* 360 degree wrap-around *)

    for w := 0 to 360 do
    begin
      tilt := slant_angle[w];

                  (** Convert tilt angle to
              an array index offset **)

      y := round(100.0 * tilt);

      (* Copy shifted el_pattern[10001] field
         values into elevation_pattern[361][1001]
         at the corresponding azimuth, downsampling
         (averaging) along the way in chunks of 10. *)

      x := y;
      z := 0;
      while (z < 1000) do
      begin
        sum := 0.0;
        for a := 0 to 9 do
        begin
          b := a + x;

          if (b >= 0) and (b <= 10000) then
            sum := sum + el_pattern[b];
          if (b < 0) then
            sum := sum + el_pattern[0];
          if (b > 10000) then
            sum := sum + el_pattern[10000];
        end;

        elevation_pattern[w][z] := sum / 10.0;
        Inc(x, 10);
        Inc(z);
      end;
    end;

    got_elevation_pattern := True;
  end;

  for x := 0 to 360 do
    for y := 0 to 1000 do
    begin
      if (got_elevation_pattern) then
        elevation := elevation_pattern[x][y]
      else
        elevation := 1.0;

      if (got_azimuth_pattern) then
        az := azimuth_pattern[x]
      else
        az := 1.0;

      FAntenna_pattern[x][y] := az * elevation;
    end;
end;


{ TPathloss }
constructor TPathloss.Create;
begin
  FLock := TCriticalSection.Create;

  FComplete := TEvent.Create(nil, True, False, '');
  FAbort := TEvent.Create(nil, True, False, '');

  FSettings := TSettings.Create;
  FMaxRange := 50;//50miles
  FModel := pmITWOM3;

  HDMode := False;
end;

destructor TPathloss.Destroy;
begin
  OnComplete := nil;

  Terminate;
  WaitFor;

  SetLength(FMask, 0);
  SetLength(FSignal, 0);
  SetLength(FAltitude, 0);
  FSettings.Free;
  FLock.Free;
  FAbort.Free;
  FComplete.Free;
  inherited;
end;

procedure TPathLoss.Terminate;
begin
  FAbort.SetEvent;
end;

procedure TPathLoss.WaitFor(const ATimeOut: cardinal = INFINITE);
begin
  FComplete.WaitFor(ATimeOut);
end;

procedure TPathLoss.SetHDMode(const AValue: boolean);
begin
  if AValue <> FHDMode then
  begin
    FHDMode := AValue;
    SetMaxRange(MaxRange);
  end;
end;

function TPathLoss.GetPixelPerDegree: integer;
begin
  if FHDMode then
    Result := 3600
  else
    Result := 1200;
end;

procedure TPathLoss.SetMaxRange(const AValue: double);
var
  size, newHeight: integer;
begin
  FLock.Enter;
  try
    FMaxRange := AValue;
    FSource.lat := 0;
    FSource.lon := 0;

    newHeight := round(2 * (FMaxRange / EARTHRADIUS_MILES) / DEG2RAD *
      PixelPerDegree);

    if newHeight <> FHeight then
    begin
      FHeight := newHeight;
      size := FHeight * FHeight;

      SetLength(FMask, size);
      SetLength(FSignal, size);
      SetLength(FAltitude, size);

      fillchar(FMask[0], size, 0);
      fillchar(FSignal[0], size, 0);
    end;
  finally
    FLock.Leave;
  end;
end;

function TPathLoss.SiteToOffset(const ALat, ALon: double;
  out X, Y, Offset: integer): boolean;
begin
  FLock.Enter;
  try
    y := FHeight - trunc(FHeight * (Alat - FMinNorth) / (FMaxNorth - FMinNorth));
    x := trunc(FHeight * ((ALon - FMinWest) / (FMaxWest - FMinWest)));

    Result := (x >= 0) and (x < FHeight) and (y >= 0) and (y < FHeight);
    if Result then
      offset := y * FHeight + x
    else
      offset := 0;
  finally
    FLock.Leave;
  end;
end;

function TPathLoss.SiteToOffset(const ALat, ALon: double; out Offset: integer): boolean;
var
  x, y: integer;
begin
  Result := SiteToOffset(ALat, ALon, x, y, Offset);
end;

function TPathloss.GetSignalValue(const AOffset: integer): byte;
begin
  Result := FSignal[AOffset];
end;

function TPathloss.GetMaskValue(const AOffset: integer): byte;
begin
  Result := FMask[AOffset];
end;

function TPathloss.PutMask(const ALat, ALon: double; const AValue: integer): integer;
var
  offset: integer;
begin
  FLock.Enter;
  try
    if SiteToOffset(ALat, ALon, offset) then
    begin
      FMask[offset] := AValue;
      Result := FMask[offset];
    end
    else
      Result := -1;
  finally
    FLock.Leave;
  end;
end;

function TPathloss.OrMask(const ALat, ALon: double; const AValue: integer): integer;
var
  offset: integer;
begin
  if not SiteToOffset(ALat, ALon, offset) then
  begin
    Result := -1;
    exit;
  end;

  FLock.Enter;
  try
    FMask[offset] := FMask[offset] or AValue;
    Result := FMask[offset];
  finally
    FLock.Leave;
  end;
end;

function TPathloss.GetMask(const ALat, ALon: double): integer;
begin
  Result := (OrMask(Alat, Alon, 0));
end;

procedure TPathloss.PutSignal(const ALat, ALon: double; const ASignal: byte);
var
  offset: integer;
begin
  (* This function writes a signal level (0-255)
     at the specified location for later recall. *)
  FLock.Enter;
  try
    if (ASignal > FHottest) then  // dBm, dBuV
      FHottest := ASignal;

    if SiteToOffset(ALat, ALon, offset) then
    begin
      // Write values to file
      FSignal[offset] := ASignal;
      exit;
    end;
  finally
    FLock.Leave;
  end;
end;

function TPathloss.GetSignal(const ALat, ALon: double): byte;
var
  offset: integer;
begin
  FLock.Enter;
  try
    if SiteToOffset(ALat, ALon, offset) then
      Result := FSignal[offset]
    else
      Result := 0;
  finally
    FLock.Leave;
  end;
end;

procedure TPathloss.PlotLOSPath(const ASource, ADestination: TSite;
  const AMaskValue: byte);
var
  block: boolean;
  x, y: integer;
  cos_xmtr_angle, cos_test_angle, test_alt: double;
  distance, rx_alt, tx_alt: double;
  Path: TPath;
  p: TPathItem;
begin
  ReadAltitudeMap(ASource); //This should be only the first time

  Path := TPath.Create(ASource, ADestination, PixelPerDegree, @GetElevation);
  try
    for y := 0 to Path.Count - 2 do
    begin
      p := TPathItem(path[y]);
      if (p.distance > FMaxRange) then
        break;
    (* Test this point only if it hasn't been already
       tested and found to be free of obstructions. *)

      if (GetMask(p.lat, p.lon) and AMaskValue) = 0 then
      begin

        distance := FEET_PER_MILE * p.distance;
        tx_alt := earthradius + ASource.alt + TPathItem(path[0]).elevation;
        rx_alt :=
          earthradius + ADestination.alt + p.elevation;

      (* Calculate the cosine of the elevation of the
         transmitter as seen at the temp rx point. *)
        if distance = 0 then
          continue;

        cos_xmtr_angle :=
          ((rx_alt * rx_alt) + (distance * distance) - (tx_alt * tx_alt)) /
          (2.0 * rx_alt * distance);

        block := False;
        for x := y downto 0 do
        begin
          distance :=
            FEET_PER_MILE * (p.distance - TPathItem(path[x]).distance);
          if distance = 0 then
            continue;
          if TPathItem(path[x]).elevation = 0 then
            test_alt :=
              earthradius + TPathItem(path[x]).elevation
          else
            test_alt :=
              earthradius + TPathItem(path[x]).elevation + FClutter;

          cos_test_angle :=
            ((rx_alt * rx_alt) + (distance * distance) - (test_alt * test_alt)) /
            (2.0 * rx_alt * distance);

        (* Compare these two angles to determine if
           an obstruction exists.  Since we're comparing
           the cosines of these angles rather than
           the angles themselves, the following "if"
           statement is reversed from what it would
           be if the actual angles were compared. *)

          if (cos_xmtr_angle >= cos_test_angle) then
          begin
            block := True;
            break;
          end;
        end;

        if not block then
          OrMask(p.lat, p.lon, AMaskValue);
      end;
    end;
  finally
    Path.Free;
  end;
end;

function TPathloss.GetAverageHeight(const AX, AY: integer; ASize: integer = 2): integer;
var
  x, y, Count, offset: integer;
begin
  Result := 0;
  Count := 0;
  for x := ax - ASize to ax + ASize do
    for y := ay - ASize to ay + ASize do
    begin
      offset := y * FHeight + x;
      if (y >= 0) and (y < FHeight) and (x >= 0) and (x < FHeight) and
        (FAltitude[offset] <> 32768) then
      begin
        Result := Result + FAltitude[offset];
        Inc(Count);
      end;
    end;
  if Count > 0 then
    Result := trunc(Result / Count);
end;

procedure TPathloss.ReadAltitudeMap(const ASource: TSite);
// Preload the map and repair the missing spots. This is kind of slow, so lets speed it up in the future
var
  dpp, lat, lon: double;
  x, y, offset: integer;
  SRTM: TSRTM;
  Elevation: double;
begin
  FLock.Enter;
  try
    if (ASource.Lat = FSource.Lat) and (ASource.Lon = FSource.Lon) then
      exit;
    FSource := ASource;

    SRTM := TSRTM.Create;
    try
      SRTM.Load(ASource.Lat, ASource.Lon, FMaxRange + RENDER_AROUND_RANGE);
      dpp := 1 / PixelPerDegree;

      y := 0;
      repeat
        lat := MinNorth + (dpp * y);
        Inc(y);
        x := 0;
        repeat
          lon := MinWest + (dpp * x);
          if SiteToOffset(lat, lon, offset) then
            FAltitude[offset] := SRTM.GetElevation(Lat, Lon);
          Inc(x);
        until lon > MaxWest;
      until lat > MaxNorth;

      for y := 0 to FHeight - 1 do
        for x := 0 to FHeight - 1 do
        begin
          offset := y * FHeight + x;
          if FAltitude[offset] = 32768 then
            FAltitude[offset] := GetAverageHeight(x, y, 2);

          Elevation := 3.28084 * FAltitude[offset];
          if Elevation > FMaxElevation then
            FMaxElevation := Elevation;
          if Elevation < FMinElevation then
            FMinElevation := Elevation;
        end;
    finally
      SRTM.Free;
    end;
  finally
    FLock.Leave;
  end;
end;

procedure TPathloss.GetElevation(Sender: TObject; ALat, ALon: double;
  out Elevation: double);
var
  offset: integer;
begin
  if SiteToOffset(ALat, ALon, offset) then
  begin
    Elevation := 3.28084 * FAltitude[offset];

    if Elevation > FMaxElevation then
      FMaxElevation := Elevation;
    if Elevation < FMinElevation then
      FMinElevation := Elevation;
  end;
end;

procedure TPathloss.PlotPropPath(const ASource, ADestination: TSite;
  const AMaskValue: byte);
(* This function plots the RF path loss between source and
          destination points based on the ITWOM propagation model,
          taking into account antenna pattern data, if available. *)
var
  x, y, ifs, ofs, errnum: integer;
  block: boolean;
  strmode: string;
  Xmtr_alt: double;
  dest_alt, xmtr_alt2, dest_alt2, cos_rcvr_angle, cos_test_angle,
  test_alt, delevation, pattern, dBm, distance, four_thirds_earth,
  rxp, field_strength: double;
  Path: TPath;
  elev: array of double;
  Temp: TSite;

  azimuth: double;
  loss: double;
  dkm: double;
  p: TPathItem;
begin
  ReadAltitudeMap(ASource); //This should be only the first time

  block := False;
  pattern := 0.0;
  cos_test_angle := 0.0;
  delevation := 0.0;
  distance := 0.0;
  field_strength := 0.0;

  Path := TPath.Create(ASource, ADestination, PixelPerDegree, @GetElevation);
  try
    setlength(elev, path.Count + 2);
    four_thirds_earth := FOUR_THIRDS * EARTHRADIUS;

    (* Copy elevations plus clutter along path into the elev[] array. *)
    for x := 1 to path.Count - 2 do
    begin
      elev[x + 2] := TPathItem(Path[x]).elevation * METERS_PER_FOOT;
      if elev[x + 2] > 0 then
        elev[x + 2] := elev[x + 2] + FClutter * METERS_PER_FOOT;
    end;

    (* Copy ending points without clutter *)

    elev[2] := TPathItem(path[0]).elevation * METERS_PER_FOOT;
    elev[Path.Count + 1] := TPathItem(path[path.Count - 1]).elevation * METERS_PER_FOOT;

  (* Since the only energy the propagation model considers
     reaching the destination is based on what is scattered
     or deflected from the first obstruction along the path,
     we first need to find the location and elevation angle
     of that first obstruction (if it exists).  This is done
     using a 4/3rds Earth radius to match the radius used by
     the irregular terrain propagation model.  This information
     is required for properly integrating the antenna's elevation
     pattern into the calculation for overall path loss. *)

    y := 2;
    while y < path.Count - 1 do
    begin
      p := TPathItem(path[y]);
      if p.distance > FMaxRange then
        break;

      if (GetMask(p.lat, p.lon) and AMaskValue) = 0 then
      begin
        distance := 5280.0 * TPathItem(path[y]).distance;
        Xmtr_alt := four_thirds_earth + ASource.alt + TPathItem(Path[0]).elevation;
        dest_alt := four_thirds_earth + ADestination.alt + TPathItem(path[y]).elevation;
        dest_alt2 := dest_alt * dest_alt;
        xmtr_alt2 := Xmtr_alt * Xmtr_alt;

      (* Calculate the cosine of the elevation of
         the receiver as seen by the transmitter. *)

        cos_rcvr_angle := ((xmtr_alt2) + (distance * distance) - (dest_alt2)) /
          (2.0 * Xmtr_alt * distance);

        if (cos_rcvr_angle > 1.0) then
          cos_rcvr_angle := 1.0;

        if (cos_rcvr_angle < -1.0) then
          cos_rcvr_angle := -1.0;

        if FGotElevationPattern then
        begin
        (* Determine the elevation angle to the first obstruction
           along the path IF elevation pattern data is available
           or an output (.ano) file has been designated. *)
          x := 2;
          block := False;
          while (x < y) and (not block) do
          begin
            distance := 5280.0 * TPathItem(path[x]).distance;

            if TPathItem(path[x]).elevation = 0.0 then
              test_alt := four_thirds_earth + TPathItem(path[x]).elevation
            else
              test_alt :=
                four_thirds_earth + TPathItem(path[x]).elevation + FClutter;

          (* Calculate the cosine of the elevation
             angle of the terrain (test point)
             as seen by the transmitter. *)

            cos_test_angle := ((xmtr_alt2) + (distance * distance) -
              (test_alt * test_alt)) / (2.0 * Xmtr_alt * distance);

            if (cos_test_angle > 1.0) then
              cos_test_angle := 1.0;

            if (cos_test_angle < -1.0) then
              cos_test_angle := -1.0;

          (* Compare these two angles to determine if
             an obstruction exists.  Since we're comparing
             the cosines of these angles rather than
             the angles themselves, the sense of the
             following "if" statement is reversed from
               what it would be if the angles themselves
             were compared. *)

            if (cos_rcvr_angle >= cos_test_angle) then
              block := True;
            Inc(x);
          end;

          if (block) then
            delevation := ((arccos(cos_test_angle)) / DEG2RAD) - 90.0
          else
            delevation := ((arccos(cos_rcvr_angle)) / DEG2RAD) - 90.0;
        end;

      (* Determine attenuation for each point along
         the path using ITWOM's point_to_point mode
         starting at y=2 (number_of_points = 1), the
         shortest distance terrain can play a role in
         path loss. *)

        elev[0] := y - 1;  (* (number of points - 1) *)

        (* Distance between elevation samples *)

        elev[1] := METERS_PER_MILE *
          (p.distance - TPathItem(path[y - 1]).distance);

        dkm := (elev[1] * elev[0]) / 1000;  // km

        case FModel of
          pmLongleyRice:
          begin
            point_to_point_ITM(elev, ASource.alt * METERS_PER_FOOT,
              ADestination.alt * METERS_PER_FOOT, FSettings.eps_dielect,
              FSettings.sgm_conductivity, FSettings.eno_ns_surfref, FSettings.frq_mhz,
              FSettings.radio_climate, FSettings.pol, FSettings.conf,
              FSettings.rel, loss,
              strmode, errnum);
          end;
          pmHATA:
          begin
            loss :=
              HATApathLoss(FSettings.frq_mhz, ASource.alt *
              METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
          end;
          pmECC33:
          begin
            loss :=
              ECC33pathLoss(FSettings.frq_mhz, ASource.alt *
              METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
          end;
          pmSUI:
          begin
            loss :=
              SUIpathLoss(FSettings.frq_mhz, ASource.alt * METERS_PER_FOOT,
              (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
          end;
          pmCost231:
          begin
            loss :=
              COST231pathLoss(FSettings.frq_mhz, ASource.alt *
              METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
          end;
          pmFreespace:
          begin
            loss := FSPLpathLoss(Settings.frq_mhz, dkm);
          end;
          pmITWOM3:
          begin
            point_to_point(elev, ASource.alt * METERS_PER_FOOT,
              ADestination.alt * METERS_PER_FOOT, FSettings.eps_dielect,
              FSettings.sgm_conductivity, FSettings.eno_ns_surfref, FSettings.frq_mhz,
              FSettings.radio_climate, FSettings.pol, FSettings.conf,
              FSettings.rel, loss,
              strmode, errnum);
          end;
          pmEricsson:
          begin
            loss :=
              EricssonpathLoss(FSettings.frq_mhz, ASource.alt *
              METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
          end;
          pmPlainEarth:
          begin
            // Plane earth
            loss := PlaneEarthLoss(dkm, ASource.alt * METERS_PER_FOOT,
              (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT));
          end;
          pmEgli:
          begin
            // Egli VHF/UHF
            loss := EgliPathLoss(FSettings.frq_mhz, ASource.alt *
              METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
              (ADestination.alt * METERS_PER_FOOT), dkm);
          end;
          pmSoil:
          begin
            // Soil
            loss := SoilPathLoss(FSettings.frq_mhz, dkm, FSettings.eps_dielect);
          end;
        end;

        temp.lat := p.lat;
        temp.lon := p.lon;

        azimuth := (GetAzimuth(ASource, temp));

      (* Integrate the antenna's radiation
         pattern into the overall path loss. *)

        x := round(10.0 * (10.0 - delevation));

        if (x >= 0) and (x <= 1000) then
        begin
          azimuth := round(azimuth);

          pattern := FSettings.antenna_pattern[trunc(azimuth)][x];

          if (pattern <> 0.0) then
          begin
            pattern := 20.0 * log10(pattern);
            loss := loss - pattern;
          end;
        end;

        if (FSettings.erp <> 0.0) then
        begin
          if (FUsedbm) then
          begin
            (* dBm is based on EIRP (ERP + 2.14) *)
            rxp := FSettings.erp / (power(10.0, (loss - 2.14) / 10.0));
            dBm := 10.0 * (log10(rxp * 1000.0));

            ifs := 200 + round(dBm);
          end
          else
          begin
            field_strength := (139.4 + (20.0 * log10(FSettings.frq_mhz)) - loss) +
              (10.0 * log10(FSettings.erp / 1000.0));

            ifs := 100 + round(field_strength);
          end;
        end
        else
        begin
          if (loss > 255) then
            ifs := 255
          else
            ifs := round(loss);
        end;

        if (ifs < 0) then
          ifs := 0;

        if (ifs > 255) then
          ifs := 255;

        ofs := GetSignal(p.lat, p.lon);
        if (ofs > ifs) then
          ifs := ofs;

        PutSignal(p.lat, p.lon, ifs);
        (* Mark this point as having been analyzed *)
        PutMask(p.lat, p.lon, int64(GetMask(p.lat, p.lon) and 7) + (AMaskValue shl 3));
      end;
      Inc(y);
    end;
  finally
    Path.Free;
  end;
end;


procedure TPathloss.Propagate(const ARange: TRangeSettings);
var
  dminwest: double;
  lon, lat: double;
  y: integer;
  dpp: double;
  edge: TSite;
begin
  FLock.Enter;
  try
    Inc(FTaskCount);
  finally
    FLock.Leave;
  end;

  dpp := 1 / PixelPerDegree;
  dminwest := dpp + ARange.min_west;
  if ARange.eastwest then
    lon := dminwest
  else
    lon := ARange.min_west;
  lat := ARange.min_north;
  y := 0;

  repeat
    if FAbort.WaitFor(0) = TWaitResult.wrSignaled then
      break;

    if (lon >= 360.0) then
      lon := lon - 360.0;

    edge.lat := lat;
    edge.lon := lon;
    edge.alt := ARange.altitude;

    if (ARange.los) then
      PlotLOSPath(ARange.Source, edge, ARange.mask_value)
    else
      PlotPropPath(ARange.Source, edge, ARange.mask_value);

    Inc(y);
    if ARange.eastwest then
      lon := dminwest + (dpp * y)
    else
      lat := ARange.min_north + (dpp * y);

  until ((ARange.eastwest) and ((LonDiff(lon, ARange.max_west) > 0.0)) or
      ((not ARange.eastwest) and (lat >= ARange.max_north)));
end;

function TPathLoss.IsCalulating: boolean;
begin
  FLock.Enter;
  try
    Result := FTaskCount > 0;
  finally
    FLock.Leave;
  end;
end;

procedure TPathloss.Calculate(const ASource: TSite; const AAltitude: double;
  const ALineOfSight: boolean = False);
var
  range_min_west, range_min_north, range_max_west, range_max_north: array[0..3] of
  double;
  range: TRangeSettings;
  i: integer;
begin
  (* This function performs a 360 degree sweep around the
     transmitter site (source location), and plots the
     line-of-sight coverage of the transmitter on the ss
     generated topographic map based on a receiver located
     at the specified altitude (in feet AGL).  Results
     are stored in memory, and written out in the form
     of a topographic map when the WritePPM() function
     is later invoked. *)
  FLock.Enter;
  try
    FMinNorth := ASource.Lat - ((FMaxRange + RENDER_AROUND_RANGE) /
      EARTHRADIUS_MILES) / DEG2RAD;
    FMaxNorth := ASource.Lat + ((FMaxRange + RENDER_AROUND_RANGE) /
      EARTHRADIUS_MILES) / DEG2RAD;
    FMinWest := ASource.Lon - ((FMaxRange + RENDER_AROUND_RANGE) / EARTHRADIUS_MILES) /
      DEG2RAD / Cos(ASource.Lat * DEG2RAD);
    FMaxWest := ASource.Lon + ((FMaxRange + RENDER_AROUND_RANGE) / EARTHRADIUS_MILES) /
      DEG2RAD / Cos(ASource.Lat * DEG2RAD);
  finally
    FLock.Leave;
  end;

  FMaskValue := 1;

  // Four sections start here
  // Process north edge east/west, east edge north/south,
  // south edge east/west, west edge north/south
  range_min_west[0] := FMinWest;
  range_min_west[1] := FMinWest;
  range_min_west[2] := FMinWest;
  range_min_west[3] := FMaxWest;

  range_min_north[0] := FMaxNorth;
  range_min_north[1] := FMinNorth;
  range_min_north[2] := FMinNorth;
  range_min_north[3] := FMinNorth;

  range_max_west[0] := FMaxWest;
  range_max_west[1] := FMinWest;
  range_max_west[2] := FMaxWest;
  range_max_west[3] := FMaxWest;

  range_max_north[0] := FMaxNorth;
  range_max_north[1] := FMaxNorth;
  range_max_north[2] := FMinNorth;
  range_max_north[3] := FMaxNorth;

  FMinElevation := 32768;
  FMaxElevation := -32768;

  FTaskCount := 0;
  FComplete.ResetEvent;
  FAbort.ResetEvent;
  for i := 0 to 3 do
  begin
    range.los := ALineOfSight;
    range.eastwest := range_min_west[i] <> range_max_west[i];
    range.min_west := range_min_west[i];
    range.max_west := range_max_west[i];
    range.min_north := range_min_north[i];
    range.max_north := range_max_north[i];

    range.altitude := AAltitude;
    range.Source := ASource;
    range.mask_value := FMaskValue;
    TProcessThread.Create(range, self);
    //Propagate(range);
  end;

  case FMaskValue of
    1: FMaskValue := 8;
    8: FMaskValue := 16;
    16: FMaskValue := 32;
    else
  end;
end;

procedure TPathLoss.CompleteTask;
var
  bComplete: boolean;
begin
  FLock.Enter;
  try
    Dec(FTaskCount);
    bComplete := FTaskCount <= 0;
  finally
    FLock.Leave;
  end;

  if bComplete then
  begin
    FComplete.SetEvent;

    if assigned(FOnComplete) then
      FOnComplete(self);
  end;
end;


{ TCalculatorItem }
constructor TCalculatorItem.Create(const ALat, ALon: double; const ASignal: byte);
begin
  FLat := ALat;
  FLon := ALon;
  FSignal := ASignal;
end;

{ TCalculator }
constructor TCalculator.Create;
begin
  inherited Create(True);
  FSRTM := TSRTM.Create;
end;

destructor TCalculator.Destroy;
begin
  FSRTM.Free;
  inherited;
end;

procedure TCalculator.GetElevation(Sender: TObject; ALat, ALon: double;
  out Elevation: double);
begin
  Elevation := 3.28084 * FSRTM.GetElevation(ALat, ALon);
end;

function TCalculator.GetPixelPerDegree: integer;
begin
  if FHDMode then
    Result := 3600
  else
    Result := 1200;
end;

{ TLOSCalculator }
procedure TLOSCalculator.Calculate(const ASource, ADestination: TSite);
var
  block: boolean;
  x, y: integer;
  cos_xmtr_angle, cos_test_angle, test_alt: double;
  distance, rx_alt, tx_alt: double;
  Path: TPath;
  p: TPathItem;
begin
  FSRTM.Load(ASource.Lat, ASource.Lon,
    GetDistance(ASource, ADestination));

  Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
  try
    for y := 0 to Path.Count - 2 do
    begin
      p := TPathItem(path[y]);
    (* Test this point only if it hasn't been already
       tested and found to be free of obstructions. *)

      distance := FEET_PER_MILE * p.distance;
      tx_alt := earthradius + ASource.alt + TPathItem(path[0]).elevation;
      rx_alt :=
        earthradius + ADestination.alt + p.elevation;

      (* Calculate the cosine of the elevation of the
         transmitter as seen at the temp rx point. *)
      if distance = 0 then
        continue;

      cos_xmtr_angle :=
        ((rx_alt * rx_alt) + (distance * distance) - (tx_alt * tx_alt)) /
        (2.0 * rx_alt * distance);

      block := False;
      for x := y downto 0 do
      begin
        distance :=
          FEET_PER_MILE * (p.distance - TPathItem(path[x]).distance);
        if distance = 0 then
          continue;
        if TPathItem(path[x]).elevation = 0 then
          test_alt :=
            earthradius + TPathItem(path[x]).elevation
        else
          test_alt :=
            earthradius + TPathItem(path[x]).elevation + FClutter;

        cos_test_angle :=
          ((rx_alt * rx_alt) + (distance * distance) - (test_alt * test_alt)) /
          (2.0 * rx_alt * distance);

        (* Compare these two angles to determine if
           an obstruction exists.  Since we're comparing
           the cosines of these angles rather than
           the angles themselves, the following "if"
           statement is reversed from what it would
           be if the actual angles were compared. *)

        if (cos_xmtr_angle >= cos_test_angle) then
        begin
          block := True;
          break;
        end;
      end;

      if not block then
        Add(TCalculatorItem.Create(p.lat, p.lon, 1));
    end;
  finally
    Path.Free;
  end;
end;


{ TPropCalculator }

constructor TPropCalculator.Create;
begin
  inherited;
  FSettings := TSettings.Create;
  FModel := TPathLossModel.pmITWOM3;
end;

destructor TPropCalculator.Destroy;
begin
  FSettings.Free;
  inherited;
end;

procedure TPropCalculator.Calculate(const ASource, ADestination: TSite);
var
  x, y, ifs, errnum: integer;
  block: boolean;
  strmode: string;
  Xmtr_alt: double;
  dest_alt, xmtr_alt2, dest_alt2, cos_rcvr_angle, cos_test_angle,
  test_alt, delevation, pattern, dBm, distance, four_thirds_earth,
  rxp, field_strength: double;
  Path: TPath;
  elev: array of double;
  Temp: TSite;

  azimuth: double;
  loss: double;
  dkm: double;
  p: TPathItem;
begin
  block := False;
  pattern := 0.0;
  cos_test_angle := 0.0;
  delevation := 0.0;
  distance := 0.0;
  field_strength := 0.0;

  FSRTM.Load(ASource.Lat, ASource.Lon,
    GetDistance(ASource, ADestination));

  Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
  try
    setlength(elev, path.Count + 2);
    four_thirds_earth := FOUR_THIRDS * EARTHRADIUS;

    (* Copy elevations plus clutter along path into the elev[] array. *)
    for x := 1 to path.Count - 2 do
    begin
      elev[x + 2] := TPathItem(Path[x]).elevation * METERS_PER_FOOT;
      if elev[x + 2] > 0 then
        elev[x + 2] := elev[x + 2] + FClutter * METERS_PER_FOOT;
    end;

    (* Copy ending points without clutter *)

    elev[2] := TPathItem(path[0]).elevation * METERS_PER_FOOT;
    elev[Path.Count + 1] := TPathItem(path[path.Count - 1]).elevation * METERS_PER_FOOT;

  (* Since the only energy the propagation model considers
     reaching the destination is based on what is scattered
     or deflected from the first obstruction along the path,
     we first need to find the location and elevation angle
     of that first obstruction (if it exists).  This is done
     using a 4/3rds Earth radius to match the radius used by
     the irregular terrain propagation model.  This information
     is required for properly integrating the antenna's elevation
     pattern into the calculation for overall path loss. *)

    y := 2;
    while y < path.Count - 1 do
    begin
      p := TPathItem(path[y]);

      distance := 5280.0 * TPathItem(path[y]).distance;
      Xmtr_alt := four_thirds_earth + ASource.alt + TPathItem(Path[0]).elevation;
      dest_alt := four_thirds_earth + ADestination.alt + TPathItem(path[y]).elevation;
      dest_alt2 := dest_alt * dest_alt;
      xmtr_alt2 := Xmtr_alt * Xmtr_alt;

      (* Calculate the cosine of the elevation of
         the receiver as seen by the transmitter. *)

      cos_rcvr_angle := ((xmtr_alt2) + (distance * distance) - (dest_alt2)) /
        (2.0 * Xmtr_alt * distance);

      if (cos_rcvr_angle > 1.0) then
        cos_rcvr_angle := 1.0;

      if (cos_rcvr_angle < -1.0) then
        cos_rcvr_angle := -1.0;

      if FGotElevationPattern then
      begin
        (* Determine the elevation angle to the first obstruction
           along the path IF elevation pattern data is available
           or an output (.ano) file has been designated. *)
        x := 2;
        block := False;
        while (x < y) and (not block) do
        begin
          distance := 5280.0 * TPathItem(path[x]).distance;

          if TPathItem(path[x]).elevation = 0.0 then
            test_alt := four_thirds_earth + TPathItem(path[x]).elevation
          else
            test_alt :=
              four_thirds_earth + TPathItem(path[x]).elevation + FClutter;

          (* Calculate the cosine of the elevation
             angle of the terrain (test point)
             as seen by the transmitter. *)

          cos_test_angle := ((xmtr_alt2) + (distance * distance) -
            (test_alt * test_alt)) / (2.0 * Xmtr_alt * distance);

          if (cos_test_angle > 1.0) then
            cos_test_angle := 1.0;

          if (cos_test_angle < -1.0) then
            cos_test_angle := -1.0;

          (* Compare these two angles to determine if
             an obstruction exists.  Since we're comparing
             the cosines of these angles rather than
             the angles themselves, the sense of the
             following "if" statement is reversed from
               what it would be if the angles themselves
             were compared. *)

          if (cos_rcvr_angle >= cos_test_angle) then
            block := True;
          Inc(x);
        end;

        if (block) then
          delevation := ((arccos(cos_test_angle)) / DEG2RAD) - 90.0
        else
          delevation := ((arccos(cos_rcvr_angle)) / DEG2RAD) - 90.0;
      end;

      (* Determine attenuation for each point along
         the path using ITWOM's point_to_point mode
         starting at y=2 (number_of_points = 1), the
         shortest distance terrain can play a role in
         path loss. *)

      elev[0] := y - 1;  (* (number of points - 1) *)

      (* Distance between elevation samples *)

      elev[1] := METERS_PER_MILE *
        (p.distance - TPathItem(path[y - 1]).distance);

      dkm := (elev[1] * elev[0]) / 1000;  // km

      case FModel of
        pmLongleyRice:
        begin
          point_to_point_ITM(elev, ASource.alt * METERS_PER_FOOT,
            ADestination.alt * METERS_PER_FOOT, FSettings.eps_dielect,
            FSettings.sgm_conductivity, FSettings.eno_ns_surfref, FSettings.frq_mhz,
            FSettings.radio_climate, FSettings.pol, FSettings.conf,
            FSettings.rel, loss,
            strmode, errnum);
        end;
        pmHATA:
        begin
          loss :=
            HATApathLoss(FSettings.frq_mhz, ASource.alt *
            METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
        end;
        pmECC33:
        begin
          loss :=
            ECC33pathLoss(FSettings.frq_mhz, ASource.alt *
            METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
        end;
        pmSUI:
        begin
          loss :=
            SUIpathLoss(FSettings.frq_mhz, ASource.alt * METERS_PER_FOOT,
            (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
        end;
        pmCost231:
        begin
          loss :=
            COST231pathLoss(FSettings.frq_mhz, ASource.alt *
            METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
        end;
        pmFreespace:
        begin
          loss := FSPLpathLoss(Settings.frq_mhz, dkm);
        end;
        pmITWOM3:
        begin
          point_to_point(elev, ASource.alt * METERS_PER_FOOT,
            ADestination.alt * METERS_PER_FOOT, FSettings.eps_dielect,
            FSettings.sgm_conductivity, FSettings.eno_ns_surfref, FSettings.frq_mhz,
            FSettings.radio_climate, FSettings.pol, FSettings.conf,
            FSettings.rel, loss,
            strmode, errnum);
        end;
        pmEricsson:
        begin
          loss :=
            EricssonpathLoss(FSettings.frq_mhz, ASource.alt *
            METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT), dkm, FSettings.Environment);
        end;
        pmPlainEarth:
        begin
          // Plane earth
          loss := PlaneEarthLoss(dkm, ASource.alt * METERS_PER_FOOT,
            (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT));
        end;
        pmEgli:
        begin
          // Egli VHF/UHF
          loss := EgliPathLoss(FSettings.frq_mhz, ASource.alt *
            METERS_PER_FOOT, (TPathItem(path[y]).elevation * METERS_PER_FOOT) +
            (ADestination.alt * METERS_PER_FOOT), dkm);
        end;
        pmSoil:
        begin
          // Soil
          loss := SoilPathLoss(FSettings.frq_mhz, dkm, FSettings.eps_dielect);
        end;
      end;

      temp.lat := p.lat;
      temp.lon := p.lon;

      azimuth := (GetAzimuth(ASource, temp));

      (* Integrate the antenna's radiation
         pattern into the overall path loss. *)

      x := round(10.0 * (10.0 - delevation));

      if (x >= 0) and (x <= 1000) then
      begin
        azimuth := round(azimuth);

        pattern := FSettings.antenna_pattern[trunc(azimuth)][x];

        if (pattern <> 0.0) then
        begin
          pattern := 20.0 * log10(pattern);
          loss := loss - pattern;
        end;
      end;

      if (FSettings.erp <> 0.0) then
      begin
        if (FUsedbm) then
        begin
          (* dBm is based on EIRP (ERP + 2.14) *)
          rxp := FSettings.erp / (power(10.0, (loss - 2.14) / 10.0));
          dBm := 10.0 * (log10(rxp * 1000.0));

          ifs := 200 + round(dBm);
        end
        else
        begin
          field_strength := (139.4 + (20.0 * log10(FSettings.frq_mhz)) - loss) +
            (10.0 * log10(FSettings.erp / 1000.0));

          ifs := 100 + round(field_strength);
        end;
      end
      else
      begin
        if (loss > 255) then
          ifs := 255
        else
          ifs := round(loss);
      end;

      if (ifs < 0) then
        ifs := 0;

      if (ifs > 255) then
        ifs := 255;

      Add(TCalculatorItem.Create(p.lat, p.lon, ifs));
      Inc(y);
    end;
  finally
    Path.Free;
  end;
end;

end.
