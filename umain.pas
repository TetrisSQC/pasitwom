unit UMain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, ExtCtrls,
  ComCtrls, BGRABitmap, upathloss;

type

  { TfrmMain }

  TfrmMain = class(TForm)
    imgRender: TImage;
    Memo1: TMemo;
    mmoMain: TMemo;
    pgMain: TPageControl;
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

{ TfrmMain }

procedure TfrmMain.FormCreate(Sender: TObject);
var
  Dest, Source: TSite;
  Calc: TPathLossCalculator;
begin
  Test;
  //54.92806814110395;//
  Source.Lat := 47.63194444;//54.92806814110395;
  Source.Lon := -122.3538889;//8.318156693302578;
  Source.Alt := 130;
  Source.Caption := 'KOMO-TV';

  Dest.Lat := Source.lat+0.4;
  Dest.Lon := Source.lon-0.1;
  Dest.Alt := 4;
  Dest.Caption := 'Home';

  Calc := TPathLossCalculator.create;
  Calc.Calculate(Source, Dest);
  memo1.lines.Assign(calc.Report);
  Calc.free;

  FPath := TPathLoss.Create;
  FPath.Model := TPathLossModel.pmITWOM3;
  FPath.UseDBm := False;
  FPath.MaxRange := 50;
  FPath.Calculate(Source, 5, false);


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

procedure TfrmMain.tmrMainTimer(Sender: TObject);
var
  Bitmap: TBGRABitmap;
  Renderer: TPathRenderer;
begin
  Bitmap := TBGRABitmap.Create;
  Renderer := TPathRenderer.Create;
  try
    if Renderer.Render(FPath, Bitmap, rtPathLoss) then //rtLineOfSight) then
    begin
      imgRender.Picture.Assign(Bitmap);
      if not FPath.IsCalulating then
      begin
        Bitmap.savetofile(ExtractFilePath(application.exename) + 'output.png');
        tmrMain.Enabled := False;
      end;
    end;
  finally
    Renderer.Free;
    Bitmap.Free;
  end;
end;

end.
