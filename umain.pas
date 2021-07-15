unit umain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, uitwom, Math;

type

  { TForm1 }

  TForm1 = class(TForm)
    Memo1: TMemo;
    procedure FormCreate(Sender: TObject);
  private

  public

  end;

var
  Form1: TForm1;

implementation

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

{ TForm1 }

procedure TForm1.FormCreate(Sender: TObject);
var
  loss: double;
  strmode: string;
  errno: integer;
begin
  point_to_point(elev, tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
    frq_mhz, radio_climate, pol, conf, rel,
    loss, strmode, errno);

  (*
point_to_point:
  elev[0] = 18.000000
  elev[1] = 70.000208
  elev[2..n] =  159.000005, 165.000005, 166.000005, 167.000005, 164.000005, 156.000005, 151.000005, 150.000005, 150.000005, 149.000005, 148.000005, 148.000005, 148.000005, 150.000005, 149.000005, 151.000005, 160.000005, 162.000005,
tht_m: 7.620000
rht_m: 6.096000
eps_dielect: 15.000000
sgm_conductivity: 0.005000
eno_ns_surfref: 301.000000
freq_mhz: 900.000000
radio climate: 5
pol: 0
conf: 0.500000
rel: 0.500000
dbloss: 98.19
strmode: 1_Hrzn_Diff
errnum: 0
*)
  memo1.Lines.add(format('dbloss: %8.2f, mode: %s, errno: %d', [loss, strmode, errno]));

  point_to_point(elev2, tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
    frq_mhz, radio_climate, pol, conf, rel,
    loss, strmode, errno);

  (* check a short path
  -----------------------------------------------------------------------
  point_to_point:
    elev[0] = 1.000000
    elev[1] = 70.000208
    elev[2..n] =  159.000005,
    elev[n..100] =  165.000005, 166.000005, 167.000005, 164.000005, 156.000005, 151.000005, 150.000005, 150.000005, 149.000005, 148.000005, 148.000005, 148.000005, 150.000005, 149.000005, 151.000005, 160.000005, 162.000005, 160.000005, 159.000005, 157.000005, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  tht_m: 7.620000
  rht_m: 6.096000
  eps_dielect: 15.000000
  sgm_conductivity: 0.005000
  eno_ns_surfref: 301.000000
  freq_mhz: 900.000000
  radio climate: 5
  pol: 0
  conf: 0.500000
  rel: 0.500000
  dbloss: 90.45
  strmode: L-o-S
  errnum: 4
  -----------------------------------------------------------------------
  *)

  memo1.Lines.add(format('dbloss: %8.2f, mode: %s, errno: %d', [loss, strmode, errno]));

  (*
   dbloss:    98,19, mode: 1_Hrzn_Diff, errno: 0
   dbloss:    90,44, mode: L-o-S, errno: 4
  *)
end;


end.
