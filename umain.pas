unit UMain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, ExtCtrls,
  ComCtrls, TAGraph, TASeries, BGRABitmap, upathloss;

type

  { TfrmMain }

  TfrmMain = class(TForm)
    crtMain: TChart;
    lsElevation: TAreaSeries;
    lsObstructions: TLineSeries;
    lsSignal: TLineSeries;
    imgRender: TImage;
    mmoReport: TMemo;
    mmoMain: TMemo;
    pgMain: TPageControl;
    tsDiagramm: TTabSheet;
    tsReport: TTabSheet;
    tsRendering: TTabSheet;
    tsTest: TTabSheet;
    tmrMain: TTimer;
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
    procedure tmrMainTimer(Sender: TObject);
  private
    FPath: TPathLoss;
  public
    procedure Test;
  end;

var
  frmMain: TfrmMain;

implementation

uses UPathRenderer, uitwom, Math;

{$R *.lfm}

const
  elev: array[0..20] of double = (18,       // num points -1
    70.000208,    // distance in meters *between elev. samples*
    // path length = 18 * 70.000208

    // a total of elev[0] + 1 heights follow:

    // [2] is the elevation at the transmitter (x=0)
    159.000005,

    // elev[0] - 1 points here:
    // [3] is the first path elev. x=elev[1] from transmitter
    165.000005, 166.000005, 167.000005, 164.000005, 156.000005,
    151.000005, 150.000005, 150.000005, 149.000005, 148.000005,
    148.000005, 148.000005, 150.000005, 149.000005, 151.000005,
    160.000005, 162.000005,

    // distance from transmitter is: x=elev[0]*elev[1]
    // [20]=[2+elev[0]] elevetion of receiver
    160.000005

    );


  tht_m = 7.62;         // transmit height
  rht_m = 6.096;         // receive height
  eps_dielect = 15;     // earth dialectric constant
  sgm_conductivity = 0.005;// earth conductivity
  eno_ns_surfref = 301;   // atmospheric refractivity at seal leavel 301 nominal
  frq_mhz = 900;       // frequency
  radio_climate = 5;       // continental temperate
  pol = 0;           // 0 = horiz, 1 = vert see itwom source
  conf = 0.50;       // 0.01 to 0.99 from itwom source
  rel = 0.50;        // 0.01 to 0.99 from itwom source


  elev2: array[0..5] of double = (1,  // number of intervals between tx site and rx site.
    70.000208,
    165,  // tx height
    // elev[0] - 1 points go here.
    165,   // rx height [elev[0]*elev[1] distance away]
    INFINITY,   // height of first interval
    NAN
    //167,   // receive height, 70 meters awawy
    //99999
    );


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

function ReadBearing(AInput: string): double;
var
  degrees, minutes, seconds: integer;
  code: integer;
begin
  Val(AInput, Result, Code);
  if Code <> 0 then
  begin
    Val(GetToken(AInput, ' '), degrees, Code);
    Val(GetToken(AInput, ' '), minutes, Code);
    Val(GetToken(AInput, ' '), seconds, Code);
    Result := Abs(degrees) + Abs(minutes / 60.0) + abs(seconds / 3600.0);
    if ((degrees < 0) or (minutes < 0) or (seconds < 0.0)) then
      Result := -Result;
    if (Result > 360.0) or (Result < -360.0) then
      Result := 0.0;
  end;
end;

(*

15.000  ; Earth Dielectric Constant (Relative permittivity)
0.005  ; Earth Conductivity (Siemens per meter)
301.000  ; Atmospheric Bending Constant (N-Units)
605.000  ; Frequency in MHz (20 MHz to 20 GHz)
5  ; Radio Climate
0  ; Polarization (0 = Horizontal, 1 = Vertical)
0.50  ; Fraction of situations
0.90  ; Fraction of time
650000  ; ERP in watts

Please consult SPLAT! documentation for the meaning and use of this data.

*)


{ TfrmMain }

procedure TfrmMain.FormCreate(Sender: TObject);
var
  Dest, Source: TSite;
  Calc: TPathLossCalculator;
  i: integer;
  az, el: TStringList;
begin
  Test;
(*
  Source.Lat := 47.63194444;//54.92806814110395;
  Source.Lon := -122.3538889;//8.318156693302578;
  Source.Alt := 130;
  Source.Caption := 'KOMO-TV';

   Source.Lat := 40.330666;
  Source.Lon := -74.120975;
  Source.Alt := 350;
  Source.Caption := 'N2SMT/R';

    //Source.Lat := 54.67463356753681;Source.Lon := 13.413096637865673;
  *)

  Source.Lon := -74.246389;//ReadBearing('40 48 8.0');
  Source.Lat := 40.802222;//ReadBearing('74 14 47.0');
  Source.Alt := 323.162750;//98.5*3.28084;
  Source.Caption := 'WNJU-DT';

  Dest.Lat := Source.lat + 0.1;
  Dest.Lon := Source.lon - 0.1;
  Dest.Alt := 4;
  Dest.Caption := 'Home';

  Calc := TPathLossCalculator.Create;
  Calc.Settings.frq_mhz := 439.250;
  Calc.Calculate(Source, Dest);
  for i := 0 to Calc.Count - 1 do
    with TCalculatorItem(Calc[i]) do
    begin
      lsSignal.AddXY(Distance, Loss);
      if IsObstruction then
        lsObstructions.AddXY(Distance, Loss);
      lsElevation.AddXY(Distance, Elevation);
    end;
  mmoReport.Lines.Assign(calc.Report);
  Calc.Free;

  FPath := TPathLoss.Create;
  FPath.Model := TPathLossModel.pmLongleyRice;//pmITWOM3;//

  az := TStringList.Create;
  el := TStringList.Create;
  az.LoadFromFile('sample_data\wnju-dt.az');
  el.LoadFromFile('sample_data\wnju-dt.el');
  FPath.Settings.LoadAntennaPattern(az, el);
  az.Free;
  el.Free;
  FPath.HDMode := False;
  FPath.UseDBm := true;
  FPath.MaxRange := 50;
  FPath.Calculate(Source, 5);
end;

procedure TfrmMain.tmrMainTimer(Sender: TObject);
var
  Bitmap: TBGRABitmap;
  Renderer: TPathRenderer;
begin
  tmrMain.Enabled := FPath.IsCalculating;
  Caption := TimeToStr(Now);
  Bitmap := TBGRABitmap.Create;
  Renderer := TPathRenderer.Create;
  try
    Renderer.UseHillShade := not tmrMain.Enabled;
    if Renderer.Render(FPath, Bitmap) then
    begin
      imgRender.Picture.Assign(Bitmap);
      if not tmrMain.Enabled then
        Bitmap.savetofile(ExtractFilePath(application.exename) + 'output.png');
    end;
  finally
    Renderer.Free;
    Bitmap.Free;
  end;
end;

procedure TfrmMain.FormDestroy(Sender: TObject);
begin
  FPath.Free;
end;

procedure TfrmMain.Test();
var
  loss: double;
  strmode: string;
  errno: integer;
begin
  mmoMain.Lines.Clear;
  (*
      dbloss: 98.19
      strmode: 1_Hrzn_Diff
      errnum: 0
  *)
  point_to_point(elev, tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
    frq_mhz, radio_climate, pol, conf, rel,
    loss, strmode, errno);

  mmoMain.Lines.add(format('dbloss: %8.2f, mode: %s, errno: %d',
    [loss, strmode, errno]));


  (*
      dbloss: 90.45
      strmode: L-o-S
      errnum: 4
  *)
  point_to_point(elev2, tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
    frq_mhz, radio_climate, pol, conf, rel,
    loss, strmode, errno);

  mmoMain.Lines.add(format('dbloss: %8.2f, mode: %s, errno: %d',
    [loss, strmode, errno]));
end;

end.
