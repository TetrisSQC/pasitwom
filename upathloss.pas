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
    FantennaPattern: TAntennaPattern;
    FHasElevationPattern: boolean;
    FHasAzimuthPattern: boolean;
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

    property HasElevationPattern: boolean read FHasElevationPattern;
    property HasAzimuthPattern: boolean read FHasAzimuthPattern;
    property AntennaPattern: TAntennaPattern read FantennaPattern;
  end;

  TPathlossModel = (pmLongleyRice, pmHATA, pmECC33, pmSUI,
    pmCost231, pmFreespace, pmITWOM3, pmEricsson, pmPlainEarth, pmEgli, pmSoil);

const
  strPropModel: array[TPathLossModel] of string = (
    'Irregular Terrain Model', 'Okumura-Hata', 'ECC33 (ITU-R P.529)',
    'Stanford University Interim', 'COST231-Hata',
    'Free space path loss (ITU-R.525)', 'ITWOM 3.0', 'Ericsson',
    'Plain Earth', 'Egli', 'Soil');

type
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
    FKnifeedge: boolean;
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
    FRenderAround: double;
    FLineOfSightOnly: Boolean;

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

    function IsCalculating: boolean;

    function SiteToOffset(const ALat, ALon: double;
      out X, Y, Offset: integer): boolean; overload;
    function SiteToOffset(const ALat, ALon: double; out Offset: integer): boolean;
      overload;

    function GetMaskValue(const AOffset: integer): byte;
    function GetSignalValue(const AOffset: integer): byte;

    procedure GetElevation(Sender: TObject; ALat, ALon: double;
      out Elevation: double); overload;
    function GetElevation(ASite: TSite): double; overload;
    procedure ReadAltitudeMap(const ASource: TSite);

    procedure Calculate(const ASource: TSite; const AAltitude: double);

    property Clutter: double read FClutter write FClutter;
    property GotElevationPattern: boolean read FGotElevationPattern
      write FGotElevationPattern;

    property Height: integer read FHeight;
    property HDMode: boolean read FHDMode write SetHDMode;
    property Model: TPathlossModel read FModel write FModel;
    property Knifeedge: boolean read FKnifeedge write FKnifeedge;

    property PixelPerDegree: integer read GetPixelPerDegree;
    property Settings: TSettings read FSettings;
    property UseDBm: boolean read FUseDBm write FUseDBM;
    property LineOfSightOnly: Boolean read FLineOfSightOnly write FLineOfSightOnly;

    property MaxRange: double read FMaxRange write SetMaxRange;
    property MinNorth: double read FMinNorth;
    property MaxNorth: double read FMaxNorth;
    property MinWest: double read FMinWest;
    property MaxWest: double read FMaxWest;
    property MaxElevation: double read FMaxElevation;
    property MinElevation: double read FMinElevation;
    property RenderAround: double read FRenderAround write FRenderAround;
    property OnComplete: TNotifyEvent read FOnComplete write FOnComplete;
  end;

  TCalculatorItem = class
  private
    FLat: double;
    FLon: double;
    FLoss: double;
    FdBm: double;
    FDistance: double;
    FIsObstruction: boolean;
    FFieldStrength: double;
    FElevation: double;
  public
    constructor Create(const ALat, ALon, ADistance, AElevation, ALoss: double;
      const ASettings: TSettings);
    property Lat: double read FLat;
    property Lon: double read FLon;
    property Loss: double read FLoss;
    property Distance: double read FDistance;
    property DBm: double read FdBM;
    property FieldStrength: double read FFieldStrength;
    property IsObstruction: boolean read FIsObstruction write FIsObstruction;
    property Elevation: Double read FElevation;
  end;

  TCalculator = class(TObjectList)
  protected
    FHDMode: boolean;
    FClutter: double;
    FSRTM: TSRTM;
    FUseMetric: boolean;
    FReport: TStrings;
    FLastElevation: double;
    procedure GetElevation(Sender: TObject; ALat, ALon: double;
      out Elevation: double); overload;
    function GetElevation(const ASite: TSite): double; overload;

    function GetPixelPerDegree: integer;
  public
    constructor Create;
    destructor Destroy; override;

    procedure Calculate(const ASource, ADestination: TSite); virtual;

    property HDMode: boolean read FHDMode write FHDMode;
    property Clutter: double read FClutter write FClutter;
    property UseMetric: boolean read FUseMetric write FUseMetric;
    property Report: TStrings read FReport;
  end;

  TLOSCalculator = class(TCalculator)
  public
    procedure Calculate(const ASource, ADestination: TSite); override;
  end;

  TPathLossCalculator = class(TCalculator)
  private
    FModel: TPathlossModel;
    FFresnelZoneClearance: double;
    FSettings: TSettings;
    FRXGain: double;

    function ElevationAngle(const ASource, ADestination: TSite): double;
    function ElevationAngle2(const ASource, ADestination: TSite): double;

    procedure SetObstruction(const ASite: TSite);
    procedure ObstructionAnalysis(const ASource, ADestination: TSite);
  public
    constructor Create;
    destructor Destroy; override;

    procedure Calculate(const ASource, ADestination: TSite); override;

    property Settings: TSettings read FSettings;
    property Model: TPathLossModel read FModel write FModel;
    property FresnelZoneClearance: double read FFresnelZoneClearance
      write FFresnelZoneClearance;
    property RXGain: double read FRXGain write FRXGain;
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
  KM_PER_MILE = 1.609344;

  EARTHRADIUS = 20902230.97;
  EARTHRADIUS_MILES = 3959.0;

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

destructor TProcessThread.Destroy;
begin
  if assigned(FParent) then
    FParent.CompleteTask();
  inherited;
end;


function PolToStr(const APolarity: integer): string;
begin
  if APolarity = 0 then
    Result := 'Horizontal'
  else
    Result := 'Vertical';
end;

function ClimateToStr(const AClimate: integer): string;
begin
  case AClimate of
    1: Result := 'Equatorial';
    2: Result := 'Continental Subtropical';
    3: Result := 'Maritime Subtropical';
    4: Result := 'Desert';
    5: Result := 'Continental Temperate';
    6: Result := 'Maritime Temperate, Over Land';
    7: Result := 'Maritime Temperate, Over Sea';
    else
      Result := 'Unknown';
  end;
end;

(*
 * Acute Angle from Rx point to an obstacle of height (opp) and
 * distance (adj)
 *)
function incidenceAngle(const opp, adj: double): double;
begin
  Result := arctan2(opp, adj) * 180 / PI;
end;

(*
 * Knife edge diffraction:
 * This is based upon a recognised formula like Huygens, but trades
 * thoroughness for increased speed which adds a proportional diffraction
 * effect to obstacles.
 *)
function ked(const Elev: array of double; const freq, rxh: double; dkm: double): double;
var
  obh, obd, rxobaoi, d: double;
  n: integer;
begin
  rxobaoi := 0;

  obh := 0;    // Obstacle height
  obd := 0;    // Obstacle distance

  dkm := dkm * 1000;  // KM to metres

  // walk along path
  for n := 2 to round(dkm / elev[1]) - 1 do
  begin
    d := (n - 2) * elev[1];  // no of points * delta = km

    //Find dip(s)
    if (elev[n] < obh) then
    begin

      // Angle from Rx point to obstacle
      rxobaoi :=
        incidenceAngle((obh - (elev[n] + rxh)), d - obd);
    end
    else
    begin
      // Line of sight or higher
      rxobaoi := 0;
    end;

    //note the highest point
    if (elev[n] > obh) then
    begin
      obh := elev[n];
      obd := d;
    end;

  end;

  if (rxobaoi >= 0) then
    Result := (rxobaoi / (300 / freq)) +
      3  // Diffraction angle divided by wavelength (m)
  else
    Result := 1;
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

  (* This function returns the azimuth (in degrees) to the
     destination as seen from the location of the source. *)

var
  dest_lat, dest_lon, src_lat, src_lon, beta, azimuth, diff, num,
  den, fraction: double;
begin
  dest_lat := ADestination.lat * DEG2RAD;
  dest_lon := ADestination.lon * DEG2RAD;

  src_lat := ASource.lat * DEG2RAD;
  src_lon := ASource.lon * DEG2RAD;

  (* Calculate Surface Distance *)

  beta :=
    arccos(sin(src_lat) * sin(dest_lat) + cos(src_lat) * cos(dest_lat) *
    cos(src_lon - dest_lon));

  (* Calculate Azimuth *)

  num := sin(dest_lat) - (sin(src_lat) * cos(beta));
  den := cos(src_lat) * sin(beta);
  fraction := num / den;

  (* Trap potential problems in acos() due to rounding *)

  if (fraction >= 1.0) then
    fraction := 1.0;

  if (fraction <= -1.0) then
    fraction := -1.0;

  (* Calculate azimuth *)

  azimuth := arccos(fraction);

  (* Reference it to True North *)

  diff := dest_lon - src_lon;

  if (diff <= -PI) then
    diff := diff + TWOPI;

  if (diff >= PI) then
    diff := diff - TWOPI;

  if (diff > 0.0) then
    azimuth := TWOPI - azimuth;

  Result := (azimuth / DEG2RAD);
end;

function arccos2(x, y: double): double;
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

procedure GetDimension(ASource, ADestination: TSite;
  out MinNorth, MaxNorth, MinWest, MaxWest: double);
begin
  if ASource.Lat < ADestination.Lat then
    MinNorth := ASource.Lat
  else
    MinNorth := ADestination.Lat;

  if ASource.Lat > ADestination.Lat then
    MaxNorth := ASource.Lat
  else
    MaxNorth := ADestination.Lat;

  if ASource.Lon < ADestination.Lon then
    MinWest := ASource.Lon
  else
    MinWest := ADestination.Lon;

  if ASource.Lon > ADestination.Lon then
    MaxWest := ASource.Lon
  else
    MaxWest := ADestination.Lon;
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
  (* This function generates a sequence of latitude and
     longitude positions between source and destination
     locations along a great circle path, and stores
     elevation and distance information for points
     along that path in the "path" structure. *)
var
  elevation: double;
  azimuth, distance, lat1, lon1, beta, den, num, lat2, lon2, total_distance,
  dx, dy, path_length, miles_per_sample, samples_per_radian: double;
  lastElevation: double;
begin
  inherited Create(True);
  Elevation := 0;
  lat1 := ASource.lat * DEG2RAD;
  lon1 := ASource.lon * DEG2RAD;
  lat2 := ADestination.lat * DEG2RAD;
  lon2 := ADestination.lon * DEG2RAD;
  samples_per_radian := APixelPerDegree * 57.295833;
  azimuth := GetAzimuth(ASource, ADestination) * DEG2RAD;

  total_distance := GetDistance(ASource, ADestination);

  if (total_distance > (30.0 / APixelPerDegree)) then
  begin
    dx := samples_per_radian * arccos(cos(lon1 - lon2));
    dy := samples_per_radian * arccos(cos(lat1 - lat2));
    path_length := sqrt((dx * dx) + (dy * dy));
    miles_per_sample := total_distance / path_length;
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
    Add(TPathItem.Create(lat1, lon1, Elevation, 0));
  end;

  distance := 0.0;
  lastElevation := 0;
  while (total_distance <> 0.0) and (distance <= total_distance) do
  begin
    distance := miles_per_sample * Count;

    beta := distance / 3959.0;
    lat2 :=
      arcsin(sin(lat1) * cos(beta) + cos(azimuth) * sin(beta) * cos(lat1));
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
    {
    while (lon2 < 0.0) do
      lon2 := lon2 + TWOPI;

    while (lon2 > TWOPI) do
      lon2 := lon2 - TWOPI;
    }
    lat2 := lat2 / DEG2RAD;
    lon2 := lon2 / DEG2RAD;

    if assigned(AGetElevation) then
      AGetElevation(Self, lat2, lon2, Elevation)
    else
      Elevation := 0;

    // fix for tile gaps in multi-tile LIDAR plots
    if (Elevation = 0) and (lastElevation > 10) then
      elevation := lastElevation;
    lastElevation := Elevation;

    Add(TPathItem.Create(lat2, lon2, Elevation, distance));
  end;

  (* Make sure exact destination point is recorded at path.length-1 *)

  if assigned(AGetElevation) then
    AGetElevation(Self, lat2, lon2, Elevation)
  else
    Elevation := 0;
  Add(TPathItem.Create(ADestination.Lat, ADestination.Lon, Elevation,
    total_distance));
end;

{ TSettings }
constructor TSettings.Create;
begin
  eps_dielect := 15.0; // Farmland
  sgm_conductivity := 0.005; // Farmland
  eno_ns_surfref := 301.0;
  frq_mhz := 605.0;
  radio_climate := 5;   // continental
  pol := 0;   // horizontal
  conf := 0.5;
  rel := 0.9;
  erp := 650000;
end;

procedure TSettings.LoadAntennaPattern(const AAzimute, AElevation: TStrings);
var
  w, i, x, y: integer;
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
  helper: double;

  procedure Extract(const ALine: string; out First, Second: double);
  var
    i, istart, iend: integer;
  begin
    istart := length(ALine);
    iend := 0;
    for i := 1 to istart do
      if (ALine[i] in [' ', #9]) then
      begin
        iStart := i;
        break;
      end;
    for i := length(ALine) downto istart do
      if (ALine[i] in [' ', #9]) then
      begin
        iEnd := i;
        break;
      end;
    Val(copy(ALine, 1, iStart - 1), First, i);
    Val(copy(ALine, iend + 1), Second, i);
  end;

begin
  FHasAzimuthPattern := False;
  FHasElevationPattern := False;

  if assigned(AAzimute) and (AAzimute.Count > 0) then
  begin
    Extract(AAzimute[0], rotation, helper);
    fillchar({%H-}azimuth, length(azimuth) * sizeof(azimuth[0]), 0);
    fillchar({%H-}read_count, length(read_count) * sizeof(read_count[0]), 0);

    for i := 1 to AAzimute.Count - 1 do
    begin
      Extract(AAzimute[i], helper, amplitude);
      x := trunc(helper);
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

    FHasAzimuthPattern := True;
  end;

  if assigned(AElevation) and (AElevation.Count > 0) then
  begin
    (* Clear azimuth pattern array *)

    Fillchar({%H-}el_pattern, length(el_pattern) * sizeof(el_pattern[0]), 0);
    Fillchar({%H-}read_count, length(read_count) * sizeof(read_count[0]), 0);

    extract(AElevation[0], mechanical_tilt, tilt_azimuth);
    for i := 0 to AElevation.Count - 1 do
    begin
      extract(AElevation[i], elevation, amplitude);
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

    FHasElevationPattern := True;
  end;

  for x := 0 to 360 do
    for y := 0 to 1000 do
    begin
      if (FHasElevationPattern) then
        elevation := elevation_pattern[x][y]
      else
        elevation := 1.0;

      if (FHasAzimuthPattern) then
        az := azimuth_pattern[x]
      else
        az := 1.0;

      FAntennaPattern[x][y] := az * elevation;
    end;
end;


{ TPathloss }
constructor TPathloss.Create;
begin
  FLock := TCriticalSection.Create;

  FLineOfSightOnly := false;
  FRenderAround := 5;
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
    // x := FHeight - trunc(FHeight * (Alat - FMinNorth) / (FMaxNorth - FMinNorth));
    // y := trunc(FHeight * ((LonDiff(FMaxWest, ALon)) / (FMaxWest - FMinWest)));

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

  (* This function analyzes the path between the source and
     destination locations.  It determines which points along
     the path have line-of-sight visibility to the source.
     Points along with path having line-of-sight visibility
     to the source at an AGL altitude equal to that of the
     destination location are stored by setting bit 1 in the
     mask[][] array, which are displayed in green when PPM
     maps are later generated by ss. *)
var
  block: boolean;
  x, y: integer;
  cos_xmtr_angle, cos_test_angle, test_alt: double;
  distance, rx_alt, tx_alt: double;
  Path: TPath;
  xp, p: TPathItem;
begin
  ReadAltitudeMap(ASource);
  Path := TPath.Create(ASource, ADestination, PixelPerDegree, @GetElevation);
  try
    for y := 0 to path.Count - 1 do
    begin
      p := TPathItem(path[y]);
      if p.distance > FMaxRange then
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

        cos_xmtr_angle :=
          ((rx_alt * rx_alt) + (distance * distance) - (tx_alt * tx_alt)) /
          (2.0 * rx_alt * distance);

        block := False;
        for x := y downto 0 do
        begin
          xp := TPathItem(path[x]);

          distance := FEET_PER_MILE * (p.distance - xp.distance);
          test_alt := earthradius + xp.Elevation;
          if xp.Elevation > 0 then
            test_alt := test_alt + FClutter;

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

        if (not block) then
          OrMask(p.lat, p.lon, AMaskValue);
      end;
    end;
  finally
    Path.Free;
  end;
end;

procedure TPathloss.ReadAltitudeMap(const ASource: TSite);
// Preload the map and repair the missing spots. This is kind of slow, so lets speed it up in the future
var
  dpp, lat, lon: double;
  x, y, offset: integer;
  SRTM: TSRTM;
  Elevation, FLastAlt: integer;
begin
  FLock.Enter;
  try
    if (ASource.Lat = FSource.Lat) and (ASource.Lon = FSource.Lon) then
      exit;
    FSource := ASource;

    FLastAlt := 0;
    SRTM := TSRTM.Create;
    try
      SRTM.Load(ASource.Lat, ASource.Lon, FMaxRange + FRenderAround);
      dpp := 1 / PixelPerDegree;

      y := 0;
      repeat
        lat := MinNorth + (dpp * y);
        Inc(y);
        x := 0;
        repeat
          lon := MinWest + (dpp * x);
          if SiteToOffset(lat, lon, offset) then
          begin
            Elevation := SRTM.GetElevation(Lat, Lon);

            if (Elevation >= 10000) then
             Elevation := FLastAlt
            else
             FLastAlt := Elevation;

            FAltitude[offset] := Elevation;
          end;
          Inc(x);
        until lon > MaxWest;
      until lat > MaxNorth;
    finally
      SRTM.Free;
    end;

  finally
    FLock.Leave;
  end;
end;

function TPathloss.GetElevation(ASite: TSite): double;
begin
  GetElevation(self, ASite.Lat, ASite.Lon, Result);
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
var
  x, y, ifs, ofs, errnum: integer;
  block: boolean;
  strmode: string;
  loss, azimuth, pattern, xmtr_alt, dest_alt, xmtr_alt2, dest_alt2,
  cos_rcvr_angle, cos_test_angle, test_alt, elevation, distance,
  four_thirds_earth, field_strength, rxp, dBm, diffloss: double;
  temp: TSite;
  dkm: double;
  path: TPath;
  e: double;
  elev: array of double;
  p: TPathitem;
begin
  ReadAltitudeMap(ASource); //This should be only the first time

  block := False;
  pattern := 0.0;
  cos_test_angle := 0.0;
  elevation := 0.0;
  distance := 0.0;
  field_strength := 0.0;

  Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
  try

    four_thirds_earth := FOUR_THIRDS * EARTHRADIUS;
    SetLength(elev, path.Count + 2);

    (* Copy elevations plus clutter along path into the elev[] array. *)
    for x := 1 to path.Count - 2 do
    begin
      elev[x + 2] := TPathItem(Path[x]).elevation * METERS_PER_FOOT;
      if elev[x + 2] > 0 then
        elev[x + 2] := elev[x + 2] + FClutter * METERS_PER_FOOT;
    end;

    (* Copy ending points without clutter *)

    elev[2] := TPathItem(path[0]).elevation * METERS_PER_FOOT;

    elev[path.Count + 1] :=
      TPathItem(path[path.Count - 1]).Elevation * METERS_PER_FOOT;

  (* Since the only energy the Longley-Rice model considers
     reaching the destination is based on what is scattered
     or deflected from the first obstruction along the path,
     we first need to find the location and elevation angle
     of that first obstruction (if it exists).  This is done
     using a 4/3rds Earth radius to match the model used by
     Longley-Rice.  This information is required for properly
     integrating the antenna's elevation pattern into the
     calculation for overall path loss. *)

    for y := 2 to path.Count - 2 do
    begin
      p := TPathItem(path[y]);
      if p.distance > FMaxRange then
        break;

        (* Process this point only if it
          has not already been processed. *)

      if (int64(GetMask(p.lat, p.lon) and 248) <> (AMaskValue shl 3)) then
      begin
        distance := FEET_PER_MILE * p.distance;
        xmtr_alt :=
          four_thirds_earth + ASource.alt + TPathItem(path[0]).elevation;
        dest_alt :=
          four_thirds_earth + Adestination.alt + p.elevation;
        dest_alt2 := dest_alt * dest_alt;
        xmtr_alt2 := xmtr_alt * xmtr_alt;

      (* Calculate the cosine of the elevation of
         the receiver as seen by the transmitter. *)

        cos_rcvr_angle :=
          ((xmtr_alt2) + (distance * distance) - (dest_alt2)) /
          (2.0 * xmtr_alt * distance);

        if (cos_rcvr_angle > 1.0) then
          cos_rcvr_angle := 1.0;

        if (cos_rcvr_angle < -1.0) then
          cos_rcvr_angle := -1.0;

        if FSettings.HasElevationPattern then
        begin
        (* Determine the elevation angle to the first obstruction
           along the path IF elevation pattern data is available
           or an output (.ano) file has been designated. *)
          block := False;
          for x := 2 to y - 1 do
            if not block then
            begin
              distance := FEET_PER_MILE * TPathItem(path[x]).distance;

              e := TPathItem(path[x]).elevation;
              test_alt := four_thirds_earth + e;
              if e > 0 then
                test_alt := test_alt + FClutter;

          (* Calculate the cosine of the elevation
             angle of the terrain (test point)
             as seen by the transmitter. *)

              cos_test_angle :=
                ((xmtr_alt2) + (distance * distance) -
                (test_alt * test_alt)) / (2.0 * xmtr_alt * distance);

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
            end;

          if (block) then
            elevation :=
              ((arccos(cos_test_angle)) / DEG2RAD) - 90.0
          else
            elevation :=
              ((arccos(cos_rcvr_angle)) / DEG2RAD) - 90.0;
        end;

      (* Determine attenuation for each point along the
         path using a prop model starting at y=2 (number_of_points = 1), the
         shortest distance terrain can play a role in
         path loss. *)

        elev[0] := y - 1;  (* (number of points - 1) *)

        (* Distance between elevation samples *)

        elev[1] :=
          METERS_PER_MILE * (p.distance - TPathItem(path[y - 1]).distance);

        if (p.elevation < 1) then
          p.elevation := 1;

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


        if (FKnifeedge) and (FModel <> pmLongleyRice) then
        begin
          diffloss :=
            ked(elev, FSettings.frq_mhz, ADestination.alt * METERS_PER_FOOT, dkm);
          loss := loss + (diffloss);  // ;)
        end;

        //Key stage. Link dB for p2p is returned as 'loss'.

        temp.lat := p.lat;
        temp.lon := p.lon;

        azimuth := (GetAzimuth(ASource, temp));

      (* Integrate the antenna's radiation
         pattern into the overall path loss. *)

        x := trunc(10.0 * (10.0 - elevation));

        if (x >= 0) and (x <= 1000) then
        begin
          azimuth := round(azimuth);
          pattern :=
            FSettings.AntennaPattern[trunc(azimuth)][x];

          if (pattern <> 0.0) then
          begin
            pattern := 20.0 * log10(pattern);
            loss := loss - pattern;
          end;
        end;

        if (FSettings.erp <> 0.0) then
        begin
          if FUseDBm then
          begin
            (* dBm is based on EIRP (ERP + 2.14) *)

            rxp :=
              FSettings.erp / (power(10.0, (loss - 2.14) / 10.0));

            dBm := 10.0 * (log10(rxp * 1000.0));

            (* Scale roughly between 0 and 255 *)

            ifs := 200 + trunc(round(dBm));
          end
          else
          begin
            field_strength :=
              (139.4 + (20.0 * log10(FSettings.frq_mhz)) - loss) +
              (10.0 * log10(FSettings.erp / 1000.0));

            ifs := 100 + round(field_strength);
          end;

          if (ifs < 0) then
            ifs := 0;

          if (ifs > 255) then
            ifs := 255;

          ofs :=
            GetSignal(p.lat, p.lon);

          if (ofs > ifs) then
            ifs := ofs;
        end
        else
        begin
          if (loss > 255) then
            ifs := 255
          else
            ifs := round(loss);

          ofs :=
            GetSignal(p.lat, p.lon);

          if (ofs < ifs) and (ofs <> 0) then
            ifs := ofs;
        end;

        PutSignal(p.lat, p.lon, ifs);
        PutMask(p.lat, p.lon, int64(GetMask(p.lat, p.lon) and 7) + (AMaskValue shl 3));
      end;
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

    sleep(1);
  until ((ARange.eastwest) and ((LonDiff(lon, ARange.max_west) > 0.0)) or
      ((not ARange.eastwest) and (lat >= ARange.max_north)));
end;

function TPathLoss.IsCalculating: boolean;
begin
  FLock.Enter;
  try
    Result := FTaskCount > 0;
  finally
    FLock.Leave;
  end;
end;

procedure TPathloss.Calculate(const ASource: TSite; const AAltitude: double);
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
    FMinNorth := ASource.Lat - ((FMaxRange + FRenderAround) /
      EARTHRADIUS_MILES) / DEG2RAD;
    FMaxNorth := ASource.Lat + ((FMaxRange + FRenderAround) /
      EARTHRADIUS_MILES) / DEG2RAD;
    FMinWest := ASource.Lon - ((FMaxRange + FRenderAround) / EARTHRADIUS_MILES) /
      DEG2RAD / Cos(ASource.Lat * DEG2RAD);
    FMaxWest := ASource.Lon + ((FMaxRange + FRenderAround) / EARTHRADIUS_MILES) /
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
    range.los := FLineOfSightOnly;
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
constructor TCalculatorItem.Create(const ALat, ALon, ADistance, AElevation, ALoss: double;
  const ASettings: TSettings);
var
  rxp: double;
begin
  FLat := ALat;
  FLon := ALon;
  FDistance := ADistance;
  FIsObstruction := False;
  FLoss := ALoss;
  FElevation := AElevation;

  if assigned(ASettings) then
  begin
    rxp := ASettings.erp / (power(10.0, (FLoss - 2.14) / 10.0));
    FdBm := 10.0 * (log10(rxp * 1000.0));
    FFieldStrength := (139.4 + (20.0 * log10(ASettings.frq_mhz)) - FLoss) +
      (10.0 * log10(ASettings.erp / 1000.0));
  end;
end;

{ TCalculator }
constructor TCalculator.Create;
begin
  inherited Create(True);
  FSRTM := TSRTM.Create;
  FReport := TStringList.Create;
end;

destructor TCalculator.Destroy;
begin
  FSRTM.Free;
  FReport.Free;
  inherited;
end;

procedure TCalculator.GetElevation(Sender: TObject; ALat, ALon: double;
  out Elevation: double);
begin
  Elevation := FSRTM.GetElevation(ALat, ALon);
  if Elevation = 32768 then
    Elevation := FLastElevation
  else
    FLastElevation := Elevation;
  Elevation := 3.28084 * Elevation;
end;

function TCalculator.GetElevation(const ASite: TSite): double;
begin
  GetElevation(self, ASite.Lat, ASite.Lon, Result);
end;

function TCalculator.GetPixelPerDegree: integer;
begin
  if FHDMode then
    Result := 3600
  else
    Result := 1200;
end;

procedure TCalculator.Calculate(const ASource, ADestination: TSite);
begin
  Clear;
  FSRTM.Load(ASource.Lat, ASource.Lon, ADestination.Lat, ADestination.Lon);
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
  inherited;
  Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
  try
    for y := 0 to Path.Count - 1 do
    begin
      p := TPathItem(path[y]);
    (* Test this point only if it hasn't been already
       tested and found to be free of obstructions. *)

      distance := FEET_PER_MILE * p.distance;
      tx_alt := EARTHRADIUS + ASource.alt + TPathItem(path[0]).elevation;
      rx_alt :=
        EARTHRADIUS + ADestination.alt + p.elevation;

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
            EARTHRADIUS + TPathItem(path[x]).elevation
        else
          test_alt :=
            EARTHRADIUS + TPathItem(path[x]).elevation + FClutter;

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
        Add(TCalculatorItem.Create(p.lat, p.lon, p.distance, p.Elevation, 1, nil));
    end;
  finally
    Path.Free;
  end;
end;


{ TPathLossCalculator }

constructor TPathLossCalculator.Create;
begin
  inherited;
  FSettings := TSettings.Create;
  FFresnelZoneClearance := 0.6;
  FModel := TPathLossModel.pmITWOM3;
end;

destructor TPathLossCalculator.Destroy;
begin
  FSettings.Free;
  inherited;
end;

function TPathLossCalculator.ElevationAngle(const ASource, ADestination: TSite): double;
var
  a, b, dx: double;
begin
  a := GetElevation(ADestination) + ADestination.alt + EARTHRADIUS;
  b := GetElevation(ASource) + ASource.alt + EARTHRADIUS;

  dx := FEET_PER_MILE * GetDistance(ASource, ADestination);

  Result := ((180.0 * (arccos(((b * b) + (dx * dx) - (a * a)) / (2.0 * b * dx))) /
    PI) - 90.0);
end;

function TPathLossCalculator.ElevationAngle2(const ASource, ADestination: TSite): double;
var
  x: integer;
  source_alt, destination_alt, cos_xmtr_angle, cos_test_angle, test_alt,
  distance, source_alt2: double;
  Path: TPath;
  p: TPathItem;

  function Arccos(x: Float): Float;
  begin
    if abs(x) = 1.0 then
      if x < 0.0 then
        arccos := Pi
      else
        arccos := 0
    else
    begin
      arccos := arctan2(sqrt((1.0 - x) * (1.0 + x)), x);
    end;
  end;

begin
  (* This function returns the angle of elevation (in degrees)
     of the destination as seen from the source location, UNLESS
     the path between the sites is obstructed, in which case, the
     elevation angle to the first obstruction is returned instead.
     "er" represents the earth radius. *)

  Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
  try

    distance := FEET_PER_MILE * GetDistance(Asource, ADestination);
    source_alt := EARTHRADIUS + Asource.alt + GetElevation(Asource);
    destination_alt := EARTHRADIUS + Adestination.alt + GetElevation(Adestination);
    source_alt2 := source_alt * source_alt;

  (* Calculate the cosine of the elevation angle of the
     destination (receiver) as seen by the source (transmitter). *)

    cos_xmtr_angle :=
      ((source_alt2) + (distance * distance) - (destination_alt * destination_alt)) /
      (2.0 * source_alt * distance);

  (* Test all points in between source and destination locations to
     see if the angle to a topographic feature generates a higher
     elevation angle than that produced by the destination.  Begin
     at the source since we're interested in identifying the FIRST
     obstruction along the path between source and destination. *)

    for x := 2 to path.Count - 1 do
    begin
      p := TPathItem(Path[x]);
      distance := FEET_PER_MILE * p.Distance;

      if p.elevation = 0 then
        test_alt := EARTHRADIUS + p.Elevation
      else
        test_alt := EARTHRADIUS + p.Elevation + FClutter;

      cos_test_angle :=
        ((source_alt2) + (distance * distance) - (test_alt * test_alt)) /
        (2.0 * source_alt * distance);

    (* Compare these two angles to determine if
       an obstruction exists.  Since we're comparing
       the cosines of these angles rather than
       the angles themselves, the sense of the
       following "if" statement is reversed from
       what it would be if the angles themselves
       were compared. *)

      if (cos_xmtr_angle >= cos_test_angle) then
      begin
        Result := ((arccos(cos_test_angle)) / DEG2RAD) - 90.0;
        exit;
      end;
    end;
    Result := ((arccos(cos_xmtr_angle)) / DEG2RAD) - 90.0;
  finally
    Path.Free;
  end;
end;

procedure TPathLossCalculator.Calculate(const ASource, ADestination: TSite);
var
  x, y, errnum: integer;
  maxloss, minloss, angle1, angle2, azimuth, pattern, patterndB: double;
  total_loss, cos_xmtr_angle, cos_test_angle, source_alt, test_alt,
  dest_alt, source_alt2, dest_alt2, distance, elevation, four_thirds_earth,
  free_space_loss, eirp, voltage, rxp, power_density, dkm: double;
  helper: string;
  dBm: double;
  Path: TPath;
  elev: array of double;
  block: boolean;
  loss: double;
  strmode: string;
  field_strength: double;
begin
  inherited;
  FReport.Clear;

  maxloss := -100000.0;
  minloss := 100000.0;
  pattern := 1.0;
  patterndB := 0.0;
  total_loss := 0.0;
  cos_test_angle := 0.0;
  free_space_loss := 0.0;
  eirp := 0.0;

  four_thirds_earth := FOUR_THIRDS * EARTHRADIUS;

  FReport.Add('Transmitter site: ' + ASource.Caption);
  FReport.Add(format('Site location: %.4f, %.4f', [ASource.lat, ASource.lon]));

  if (FUseMetric) then
  begin
    FReport.Add(format('Ground elevation: %.2f meters AMSL ',
      [METERS_PER_FOOT * GetElevation(ASource)]));
    FReport.Add(format('Antenna height: %.2f meters AGL / %.2f meters AMSL',
      [METERS_PER_FOOT * Asource.alt, METERS_PER_FOOT *
      (ASource.alt + GetElevation(ASource))]));
  end
  else
  begin
    FReport.Add('Ground elevation: %.2f feet AMSL',
      [GetElevation(ASource)]);
    FReport.Add('Antenna height: %.2f feet AGL / %.2f feet AMSL',
      [ASource.alt, ASource.alt + GetElevation(ASource)]);
  end;

  azimuth := GetAzimuth(ASource, ADestination);
  angle1 := ElevationAngle(ASource, ADestination);
  angle2 := ElevationAngle2(Asource, ADestination);

  patterndB := 0;
  if (FSettings.HasElevationPattern) or (FSettings.HasAzimuthPattern) then
  begin
    x := round(10.0 * (10.0 - angle2));

    if (x >= 0) and (x <= 1000) then
      pattern :=
        FSettings.AntennaPattern[round(azimuth)][x];

    patterndB := 20.0 * log10(pattern);
  end;

  if (FUseMetric) then
    FReport.Add(format('Distance to %s: %.2f kilometers',
      [ADestination.Caption, KM_PER_MILE * GetDistance(ASource, ADestination)]))

  else
    FReport.Add(format('Distance to %s: %.2f miles',
      [ADestination.Caption, GetDistance(ASource, ADestination)]));

  FReport.Add(format('Azimuth to %s: %.2f degrees grid',
    [ADestination.Caption, Azimuth]));


  FReport.Add(format('Downtilt angle to %s: %+.4f degrees',
    [ADestination.Caption, angle1]));



  (* Receiver *)

  FReport.Add(sLineBreak + 'Receiver site: ' + ADestination.Caption + sLinebreak);
  FReport.Add(format('Site location: %.4f, %.4f', [ADestination.lat, ADestination.lon]));

  if (FUseMetric) then
  begin
    FReport.Add(format('Ground elevation: %.2f meters AMSL',
      [METERS_PER_FOOT * GetElevation(ADestination)]));
    FReport.Add(format('Antenna height: %.2f meters AGL / %.2f meters AMSL',
      [METERS_PER_FOOT * ADestination.alt, METERS_PER_FOOT *
      (ADestination.alt + GetElevation(ADestination))]));
  end
  else
  begin
    FReport.Add(format('Ground elevation: %.2f feet AMSL',
      [GetElevation(ADestination)]));
    FReport.Add(format('Antenna height: %.2f feet AGL / %.2f feet AMSL',
      [ADestination.alt, ADestination.alt + GetElevation(ADestination)]));
  end;

  if (FUseMetric) then
    FReport.Add(format('Distance to %s: %.2f kilometers',
      [ASource.Caption, KM_PER_MILE * GetDistance(Asource, Adestination)]))
  else
    FReport.Add(format('Distance to %s: %.2f miles',
      [ASource.Caption, GetDistance(ASource, ADestination)]));

  azimuth := GetAzimuth(ADestination, ASource);

  angle1 := ElevationAngle(ADestination, ASource);
  angle2 := ElevationAngle2(ADestination, ASource);

  FReport.Add(format('Azimuth to %s: %.2f degrees grid', [ASource.Caption, azimuth]));


  FReport.Add(format('Downtilt angle to %s: %.4f degrees', [ASource.Caption, angle1]));

  if (FSettings.frq_mhz > 0.0) then
  begin
    FReport.Add(sLineBreak + sLineBreak + 'Propagation model: ' +
      strPropModel[FModel]);

    FReport.Add('Model sub-type: ' + strEnvironment[FSettings.Environment]);

    FReport.Add(format('Earth''s Dielectric Constant: %.3f',
      [FSettings.eps_dielect]));
    FReport.Add(format('Earth''s Conductivity: %.3f Siemens/meter',
      [FSettings.sgm_conductivity]));
    FReport.Add(format('Atmospheric Bending Constant (N-units): %.3f ppm',
      [FSettings.eno_ns_surfref]));
    FReport.Add(format('Frequency: %.3f MHz', [FSettings.frq_mhz]));
    FReport.Add(format('Radio Climate: %d (%s)',
      [FSettings.radio_climate, ClimateToStr(FSettings.radio_climate)]));
    FReport.Add(format('Polarisation: %d (%s)', [FSettings.pol,
      PolToStr(FSettings.pol)]));

    FReport.Add(format('Fraction of Situations: %.1f%%', [FSettings.conf * 100.0]));
    FReport.Add(format('Fraction of Time: %.1f%%', [FSettings.rel * 100.0]));

    if (FSettings.erp <> 0.0) then
    begin
      FReport.Add(sLinebreak + format('Receiver gain: %.1f dBd / %.1f dBi',
        [FRXGain, FRXGain + 2.14]));
      FReport.Add('Transmitter ERP plus Receiver gain: ');

      helper := '';
      if (FSettings.erp < 1.0) then
        helper := helper + format('%.1f milliwatts', [1000.0 * FSettings.erp]);

      if (FSettings.erp >= 1.0) and (FSettings.erp < 10.0) then
        helper := helper + format('%.1f Watts', [FSettings.erp]);

      if (FSettings.erp >= 10.0) and (FSettings.erp < 10.0e3) then
        helper := helper + format('%.0f Watts', [FSettings.erp]);

      if (FSettings.erp >= 10.0e3) then
        helper := helper + format('%.3f kilowatts', [FSettings.erp / 1.0e3]);

      dBm := 10.0 * (log10(FSettings.erp * 1000.0));
      FReport.Add(format('%s (%+.2f dBm)', [helper, dBm]));
      FReport.Add(format('Transmitter ERP minus Receiver gain: %.2f dBm',
        [dBm - FRXGain]));

      (* EIRP = ERP + 2.14 dB *)

      FReport.Add('Transmitter EIRP plus Receiver gain: ');

      eirp := FSettings.erp * 1.636816521;

      helper := '';
      if (eirp < 1.0) then
        helper := helper + format('%.1f milliwatts', [1000.0 * eirp]);

      if (eirp >= 1.0) and (eirp < 10.0) then
        helper := helper + format('%.1f Watts', [eirp]);

      if (eirp >= 10.0) and (eirp < 10.0e3) then
        helper := helper + format('%.0f Watts', [eirp]);

      if (eirp >= 10.0e3) then
        helper := helper + format('%.3f kilowatts', [eirp / 1.0e3]);

      dBm := 10.0 * (log10(eirp * 1000.0));
      FReport.Add(format('%s (%+.2f dBm)', [helper, dBm]));

      // Rx gain
      FReport.Add(format('Transmitter EIRP minus Receiver gain: %.2f dBm',
        [dBm - FRXGain]));
    end;

    FReport.Add(format('Summary for the link between %s and %s:' +
      sLineBreak, [ASource.Caption, ADestination.Caption]));

    if (patterndB <> 0.0) then
      FReport.Add(format('%s antenna pattern towards %s: %.3f (%.2f dB)',
        [ASource.Caption, ADestination.Caption, pattern, patterndB]));

    Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
    (* source=TX, destination=RX *)
    try
      SetLength(elev, Path.Count + 2);
      for x := 1 to Path.Count - 1 do
      begin
        if TPathItem(Path[x]).Elevation = 0 then
          elev[x + 2] := METERS_PER_FOOT * TPathItem(Path[x]).Elevation
        else
          elev[x + 2] :=
            METERS_PER_FOOT * (FClutter + TPathItem(Path[x]).Elevation);
      end;

      (* Copy ending points without clutter *)

      elev[2] := TPathItem(path[0]).elevation * METERS_PER_FOOT;
      elev[Path.Count + 1] :=
        TPathItem(path[path.Count - 1]).elevation * METERS_PER_FOOT;

      azimuth := round(GetAzimuth(ASource, ADestination));

      for y := 2 to path.Count - 1 do
      begin (* path.length-1 avoids LR error *)
        distance := FEET_PER_MILE * TPathItem(path[y]).distance;

        source_alt := four_thirds_earth + ASource.alt + TPathItem(path[0]).Elevation;
        dest_alt := four_thirds_earth + ADestination.alt +
          TPathItem(path[y]).Elevation;
        dest_alt2 := dest_alt * dest_alt;
        source_alt2 := source_alt * source_alt;

      (* Calculate the cosine of the elevation of
         the receiver as seen by the transmitter. *)

        cos_xmtr_angle :=
          ((source_alt2) + (distance * distance) - (dest_alt2)) /
          (2.0 * source_alt * distance);

        if (FSettings.HasElevationPattern) then
        begin
        (* If an antenna elevation pattern is available, the
           following code determines the elevation angle to
           the first obstruction along the path. *)

          block := False;
          for x := 2 to y - 1 do
          begin
            distance :=
              FEET_PER_MILE * (TPathItem(path[y]).distance -
              TPathItem(path[x]).distance);
            test_alt :=
              four_thirds_earth + TPathItem(path[x]).elevation;

          (* Calculate the cosine of the elevation
             angle of the terrain (test point)
             as seen by the transmitter. *)

            cos_test_angle :=
              ((source_alt2) + (distance * distance) -
              (test_alt * test_alt)) / (2.0 * source_alt * distance);

          (* Compare these two angles to determine if
             an obstruction exists.  Since we're comparing
             the cosines of these angles rather than
             the angles themselves, the sense of the
             following "if" statement is reversed from
             what it would be if the angles themselves
             were compared. *)

            if (cos_xmtr_angle >= cos_test_angle) then
            begin
              block := True;
              break;
            end;
          end;

        (* At this point, we have the elevation angle
           to the first obstruction (if it exists). *)
        end;

      (* Determine path loss for each point along the
         path using Longley-Rice's point_to_point mode
         starting at x=2 (number_of_points = 1), the
         shortest distance terrain can play a role in
         path loss. *)

        elev[0] := y - 1;  (* (number of points - 1) *)

        (* Distance between elevation samples *)

        elev[1] :=
          METERS_PER_MILE * (TPathItem(path[y]).distance -
          TPathItem(path[y - 1]).distance);

      (*
         point_to_point(elev, source.alt*METERS_PER_FOOT,
         destination.alt*METERS_PER_FOOT, LR.eps_dielect,
         LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
         LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
         strmode, errnum);
       *)
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

        if (block) then
          elevation :=
            ((arccos(cos_test_angle)) / DEG2RAD) - 90.0
        else
          elevation :=
            ((arccos(cos_xmtr_angle)) / DEG2RAD) - 90.0;

      (* Integrate the antenna's radiation
         pattern into the overall path loss. *)

        patterndB := 0.0;
        if FSettings.HasElevationPattern then
        begin
          x := round(10.0 * (10.0 - elevation));

          if (x >= 0) and (x <= 1000) then
          begin
            pattern :=
              FSettings.AntennaPattern[round(azimuth)][x];

            if (pattern <> 0.0) then
            begin
              patterndB := 20.0 * log10(pattern);
            end;
          end;

        end;

        total_loss := loss - patterndB;

        if (total_loss > maxloss) then
          maxloss := total_loss;

        if (total_loss < minloss) then
          minloss := total_loss;

        with TPathItem(path[y]) do
        Add(TCalculatorItem.Create(lat, lon,
          Distance, Elevation, total_loss, FSettings));
      end;

      distance := GetDistance(ASource, ADestination);

      if (distance <> 0.0) then
      begin
        free_space_loss :=
          36.6 + (20.0 * log10(FSettings.frq_mhz)) + (20.0 * log10(distance));
        FReport.Add(format('Free space path loss: %.2f dB', [free_space_loss]));
      end;

      FReport.Add(format('Computed path loss: %.2f dB', [loss]));


      if ((loss * 1.5) < free_space_loss) then
      begin
        FReport.Add(format('Model error! Computed loss of %.1fdB is ' +
          'greater than free space loss of %.1fdB. ' +
          'Check your inuts for model %s', [loss, free_space_loss,
          strPropModel[FModel]]));
        exit;
      end;

      if (free_space_loss <> 0.0) then
        FReport.Add(format('Attenuation due to terrain shielding: %.2f dB',
          [loss - free_space_loss]));

      if (patterndB <> 0.0) then
      begin
        FReport.Add(format('Total path loss including %s antenna pattern: %.2f dB',
          [ASource.Caption, total_loss]));
      end;

      if (FSettings.erp <> 0.0) then
      begin
        field_strength :=
          (139.4 + (20.0 * log10(FSettings.frq_mhz)) - total_loss) +
          (10.0 * log10(FSettings.erp / 1000.0));

        (* dBm is referenced to EIRP *)

        rxp := eirp / (power(10.0, (total_loss / 10.0)));
        dBm := 10.0 * (log10(rxp * 1000.0));
        power_density :=
          (eirp / (power(10.0, (total_loss - free_space_loss) / 10.0)));
        (* divide by 4*PI*distance_in_meters squared *)
        power_density := power_density / (4.0 * PI * distance * distance * 2589988.11);

        FReport.Add(format('Field strength at %s: %.2f dBuV/meter',
          [ADestination.Caption, field_strength]));
        FReport.Add(format('Signal power level at %s: %+.2f dBm',
          [ADestination.Caption, dBm]));
        FReport.Add(format('Signal power density at %s: %+.2f dBW per square meter',
          [ADestination.Caption, 10.0 * log10(power_density)]));
        voltage :=
          1.0e6 * sqrt(50.0 * (eirp /
          (power(10.0, (total_loss - 2.14) / 10.0))));
        FReport.Add(format('Voltage across 50 ohm dipole at %s: %.2f uV (%.2f dBuV)\',
          [Adestination.Caption, voltage, 20.0 * log10(voltage)]));

        voltage :=
          1.0e6 * sqrt(75.0 * (eirp /
          (power(10.0, (total_loss - 2.14) / 10.0))));
        FReport.Add(format('Voltage across 75 ohm dipole at %s: %.2f uV (%.2f dBuV)',
          [ADestination.Caption, voltage, 20.0 * log10(voltage)]));
      end;

      if (FModel in [pmLongleyRice, pmITWOM3]) then
      begin
        case errnum of
          0: helper := ' (No error)';
          1: helper := sLinebreak +
              '  Warning: Some parameters are nearly out of range.' +
              sLinebreak + '  Results should be used with caution.';
          2: helper := sLinebreak +
              '  Note: Default parameters have been substituted for impossible ones.';
          3: helper := sLinebreak +
              '  Warning: A combination of parameters is out of range for this model.'
              + sLinebreak + '  Results should be used with caution.';
          else
            helper :=
              sLineBreak + '  Warning: Some parameters are out of range for this model.'
              + sLineBreak + '  Results should be used with caution.';
        end;
        FReport.Add(format('%s model error number: %d%s',
          [strPropModel[FModel], errnum, helper]));
      end;
    finally
      Path.Free;
    end;
  end;
  ObstructionAnalysis(ASource, Adestination);
end;


procedure TPathLossCalculator.SetObstruction(const ASite: TSite);
var
  i: integer;
begin
  for i := 0 to Count - 1 do
    with TCalculatorItem(items[i]) do
      if (Lat = ASite.Lat) and (Lon = ASite.Lon) then
      begin
        IsObstruction := True;
        exit;
      end;
end;

procedure TPathLossCalculator.ObstructionAnalysis(const ASource, ADestination: TSite);
var
  x: integer;
  site_x: TSite;
  h_r, h_t, h_x, h_r_orig, cos_tx_angle, cos_test_angle, cos_tx_angle_f1,
  cos_tx_angle_fpt6, d_tx, d_x, h_r_f1, h_r_fpt6, h_f, h_los, lambda: double;
  helper, fpt6, f1: string;
  Path: TPath;
begin
  (* Perform an obstruction analysis along the
     path between receiver and transmitter. *)

  lambda := 0.0;

  Path := TPath.Create(ASource, ADestination, GetPixelPerDegree, @GetElevation);
  try
    h_r := GetElevation(ADestination) + ADestination.alt + EARTHRADIUS;
    h_r_f1 := h_r;
    h_r_fpt6 := h_r;
    h_r_orig := h_r;
    h_t := GetElevation(ASource) + ASource.alt + EARTHRADIUS;
    d_tx := FEET_PER_MILE * GetDistance(ADestination, ASource);
    cos_tx_angle :=
      ((h_r * h_r) + (d_tx * d_tx) - (h_t * h_t)) / (2.0 * h_r * d_tx);
    cos_tx_angle_f1 := cos_tx_angle;
    cos_tx_angle_fpt6 := cos_tx_angle;

    if (FSettings.frq_mhz > 0) then
      lambda := 9.8425e8 / (FSettings.frq_mhz * 1e6);

    if (clutter > 0.0) then
    begin
      if FUseMetric then
        helper := format('%.2f meters', [METERS_PER_FOOT * clutter])
      else
        helper := format('%.2f feet', [clutter]);

      FReport.Add(format('Terrain has been raised by %s ' +
        'to account for ground clutter.', [helper]));
    end;

  (* At each point along the path calculate the cosine
     of a sort of "inverse elevation angle" at the receiver.
     From the antenna, 0 deg. looks at the ground, and 90 deg.
     is parallel to the ground.
     Start at the receiver.  If this is the lowest antenna,
     then terrain obstructions will be nearest to it.  (Plus,
     that's the way ppa!'s original los() did it.)
     Calculate cosines only.  That's sufficient to compare
     angles and it saves the extra computational burden of
     acos().  However, note the inverted comparison: if
     acos(A) > acos(B), then B > A. *)

    for x := path.Count - 1 downto 1 do
    begin
      site_x.lat := TPathItem(path[x]).lat;
      site_x.lon := TPathItem(path[x]).lon;
      site_x.alt := 0.0;

      h_x := GetElevation(site_x) + EARTHRADIUS + clutter;
      d_x := FEET_PER_MILE * GetDistance(ADestination, site_x);

      (* Deal with the LOS path first. *)

      cos_test_angle :=
        ((h_r * h_r) + (d_x * d_x) - (h_x * h_x)) / (2.0 * h_r * d_x);

      if (cos_tx_angle > cos_test_angle) then
      begin
        if (h_r = h_r_orig) then
          FReport.Add(format('Between %s and %s, ' +
            'obstructions were detected at:' + sLinebreak + sLinebreak,
            [ADestination.Caption, ASource.Caption]));

        SetObstruction(site_x);
        if (site_x.lat >= 0.0) then
        begin
          if (FUseMetric) then
            FReport.Add(format('   %8.4f N,%9.4f W, %5.2f kilometers, %6.2f meters AMSL',
              [site_x.lat, site_x.lon, KM_PER_MILE * (d_x / FEET_PER_MILE),
              METERS_PER_FOOT * (h_x - EARTHRADIUS)]))
          else
            FReport.Add(format('   %8.4f N,%9.4f W, %5.2f miles, %6.2f feet AMSL',
              [site_x.lat, site_x.lon, d_x / FEET_PER_MILE, h_x - EARTHRADIUS]));
        end
        else
        begin
          if (FUseMetric) then
            FReport.Add(format('   %8.4f S,%9.4f W, %5.2f kilometers, %6.2f meters AMSL',
              [-site_x.lat, site_x.lon, KM_PER_MILE * (d_x / FEET_PER_MILE),
              METERS_PER_FOOT * (h_x - EARTHRADIUS)]))
          else
            FReport.Add(format('   %8.4f S,%9.4f W, %5.2f miles, %6.2f feet AMSL',
              [-site_x.lat, site_x.lon, d_x / FEET_PER_MILE, h_x - EARTHRADIUS]));
        end;
      end;

      while (cos_tx_angle > cos_test_angle) do
      begin
        h_r := h_r + 1;
        cos_test_angle :=
          ((h_r * h_r) + (d_x * d_x) - (h_x * h_x)) / (2.0 * h_r * d_x);
        cos_tx_angle :=
          ((h_r * h_r) + (d_tx * d_tx) - (h_t * h_t)) / (2.0 * h_r * d_tx);
      end;

      if (FSettings.frq_mhz > 0) then
      begin
        (* Now clear the first Fresnel zone... *)

        cos_tx_angle_f1 :=
          ((h_r_f1 * h_r_f1) + (d_tx * d_tx) - (h_t * h_t)) / (2.0 * h_r_f1 * d_tx);
        h_los :=
          sqrt(h_r_f1 * h_r_f1 + d_x * d_x - 2 * h_r_f1 * d_x * cos_tx_angle_f1);
        h_f := h_los - sqrt(lambda * d_x * (d_tx - d_x) / d_tx);

        while (h_f < h_x) do
        begin
          h_r_f1 := h_r_f1 + 1;
          cos_tx_angle_f1 :=
            ((h_r_f1 * h_r_f1) + (d_tx * d_tx) - (h_t * h_t)) /
            (2.0 * h_r_f1 * d_tx);
          h_los :=
            sqrt(h_r_f1 * h_r_f1 + d_x * d_x - 2 * h_r_f1 * d_x * cos_tx_angle_f1);
          h_f :=
            h_los - sqrt(lambda * d_x * (d_tx - d_x) / d_tx);
        end;

        (* and clear the 60% F1 zone. *)

        cos_tx_angle_fpt6 :=
          ((h_r_fpt6 * h_r_fpt6) + (d_tx * d_tx) - (h_t * h_t)) /
          (2.0 * h_r_fpt6 * d_tx);
        h_los :=
          sqrt(h_r_fpt6 * h_r_fpt6 + d_x * d_x - 2 * h_r_fpt6 * d_x *
          cos_tx_angle_fpt6);
        h_f :=
          h_los - FFresnelZoneClearance * sqrt(lambda * d_x *
          (d_tx - d_x) / d_tx);

        while (h_f < h_x) do
        begin
          h_r_fpt6 := h_r_fpt6 + 1;
          cos_tx_angle_fpt6 :=
            ((h_r_fpt6 * h_r_fpt6) + (d_tx * d_tx) - (h_t * h_t)) /
            (2.0 * h_r_fpt6 * d_tx);
          h_los :=
            sqrt(h_r_fpt6 * h_r_fpt6 + d_x * d_x - 2 * h_r_fpt6 *
            d_x * cos_tx_angle_fpt6);
          h_f :=
            h_los - FFresnelZoneClearance * sqrt(lambda * d_x *
            (d_tx - d_x) / d_tx);
        end;
      end;
    end;

    if (h_r > h_r_orig) then
    begin
      if FUseMetric then
        helper := format(sLinebreak + 'Antenna at %s must be raised to ' +
          'at least %.2f meters AGL' + sLinebreak +
          'to clear all obstructions detected.',
          [ADestination.Caption, METERS_PER_FOOT *
          (h_r - GetElevation(ADestination) - EARTHRADIUS)])

      else
        helper := format(sLinebreak + 'Antenna at %s must be raised to ' +
          'at least %.2f feet AGL' + sLinebreak +
          'to clear all obstructions detected.',
          [ADestination.Caption, h_r - GetElevation(ADestination) - EARTHRADIUS]);
    end
    else
      helper := 'No obstructions to LOS path due to terrain were detected';

    if (FSettings.frq_mhz > 0) then
    begin
      if (h_r_fpt6 > h_r_orig) then
      begin
        if FUseMetric then
          fpt6 :=
            format(sLinebreak + 'Antenna at %s must be raised to at least %.2f ' +
            'meters AGL ' + sLinebreak +
            'to clear %.0f%% of the first Fresnel zone.',
            [ADestination.Caption, METERS_PER_FOOT *
            (h_r_fpt6 - GetElevation(ADestination) - EARTHRADIUS),
            FFresnelZoneClearance * 100.0])

        else
          fpt6 := format(sLinebreak + 'Antenna at %s must be raised to at least %.2f ' +
            'feet AGL' + sLinebreak +
            'to clear %.0f%% of the first Fresnel zone.',
            [ADestination.Caption, h_r_fpt6 - GetElevation(ADestination) -
            EARTHRADIUS, FFresnelZoneClearance * 100.0]);
      end
      else
        fpt6 := format(sLinebreak + '%.0f%% of the first Fresnel zone is clear.',
          [FFresnelZoneClearance * 100.0]);

      if (h_r_f1 > h_r_orig) then
      begin
        if FUseMetric then
          f1 := format(sLinebreak +
            'Antenna at %s must be raised to at least %.2f meters AGL' +
            slinebreak + 'to clear the first Fresnel zone.',
            [ADestination.Caption, METERS_PER_FOOT *
            (h_r_f1 - GetElevation(ADestination) - EARTHRADIUS)])

        else
          f1 := format(sLinebreak +
            'Antenna at %s must be raised to at least %.2f feet AGL' +
            sLinebreak + 'to clear the first Fresnel zone.',
            [ADestination.Caption, h_r_f1 - GetElevation(ADestination) - EARTHRADIUS]);

      end
      else
        f1 := slinebreak + 'The first Fresnel zone is clear.';
    end;

    FReport.add(helper);

    if (FSettings.frq_mhz > 0) then
    begin
      FReport.add(f1);
      FReport.add(fpt6);
    end;
  finally
    Path.Free;
  end;
end;

end.
