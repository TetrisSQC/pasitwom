(******************************************************************************
* ITWOM version 3.0a, January 20, 2011  File: itwom3.0a.cpp                   *
* Provenance:  Further test version of itwom2.0m re adj to Hrzn rangefactors  *
* 1. This file is based on a thorough debugging, completion, and update of    *
* the ITM, based on an original, public domain version of this file obtained  *
* from: ftp://flattop.its.bldrdoc.gov/itm/ITMDLL.cpp prior to May, 2007. C++  *
* routines for this program are taken from a translation of the FORTRAN code  *
* written by U.S. Department of Commerce NTIA/ITS Institute for               *
* Telecommunication Sciences Irregular Terrain Model (ITM) (Longley-Rice).    *
* 2. The Linux version of this file incorporates improvements suggested by a  *
* study of changes made to file itm.cpp by J. D. McDonald to remove Microsoft *
* Windows dll-isms and to debug an ambguity in overloaded calls.              *
* 3. The Linux version of this file also incorporates improvements suggested  *
* by a study of further modifications made to itm.cpp by John A. Magliacane   *
* to remove unused variables, unneeded #includes, and to replace power()      *
* statements with explicit multiplications to improve execution speed and     *
* accuracy.                                                                   *
* 4. On August 19, 2007 this file was modified by Sid Shumate to include      *
* changes and updates included in version 7.0 of ITMDLL.cpp, which was        *
* released by the NTIA/ITS on June 26, 2007. With correction set SS1 and      *
* SS2: itm71.cpp.                                                             *
* 5. On Feb. 5, 2008 this file became v.1.0 of the ITWOM with the addition,   *
* by Sid Shumate, of multiple corrections, the replacement of subroutines     *
* lrprop and alos with lrprop2 and alos2, and the addition of subroutine      *
* saalos to incorporate Radiative Transfer Engine (RTE) computations in the   *
* line of sight range.                                                        *
* Update 8 Jun 2010 to modify alos to match 2010 series of IEEE-BTS           *
* newsletter articles                                                         *
* Update June 12, 2010 to z version to change test outputs                    *
* Offshoot start date June 23, 2010 to start itwom2.0 dual version for FCC.   *
* Update to 2.0b July 25 to correct if statement errors in adiff2 re two peak *
* calculations starting at line 525                                           *
* Development to 2.0c 8 Aug 2010 after modifying saalos and adiff for full    *
* addition of saalos treatment to post obstruction calculations and           *
* debugging.                                                                  *
* Modified to make 1st obs loss=5.8 only, no clutter loss considered          *
*                                                                             *
* Commented out unused variables and calculations to eliminate gcc warnings   *
*    (-Wunused-but-set-variable)  -- John A. Magliacane -- July 25, 2013      *
*                                                                             *
* The additional models are made for Signal Server by Alex Farrant            *
* Translation do Pascal made by Christian Hackbart                            *
*******************************************************************************)

unit uitwom;

{$mode objfpc}{$H+}

interface

procedure point_to_point(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref, frq_mhz: double;
  const radio_climate, pol: integer; const conf, rel: double;
  out dbloss: double; out strmode: string; out errnum: integer);

procedure point_to_point_ITM(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref, frq_mhz: double;
  const radio_climate, pol: integer; const conf, rel: double;
  out dbloss: double; out strmode: string; out errnum: integer);

procedure point_to_pointMDH_two(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
  enc_ncc_clcref, clutter_height, clutter_density, delta_h_diff, frq_mhz: double;
  const radio_climate, pol: integer; const timepct, locpct, confpct: double;
  out dbloss: double; out propmode: integer; out deltaH: double; out errnum: integer);

procedure point_to_pointDH(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
  enc_ncc_clcref, clutter_height, clutter_density, delta_h_diff, frq_mhz: double;
  const radio_climate, pol: integer; const conf, rel: double;
  out dbloss, deltaH: double; out errnum: integer);

procedure area(const ModVar: longint; const deltaH, tht_m, rht_m, dist_km: double;
  const TSiteCriteria, RSiteCriteria: integer;
  const eps_dielect, sgm_conductivity, eno_ns_surfref, enc_ncc_clcref,
  clutter_height, clutter_density, delta_h_diff, frq_mhz: double;
  const radio_climate, pol: integer; const pctTime, pctLoc, pctConf: double;
  out dbloss: double; out errnum: integer);

function ITMAreadBLoss(const ModVar: integer;
  const deltaH, tht_m, rht_m, dist_km: double;
  const TSiteCriteria, RSiteCriteria: integer;
  const eps_dielect, sgm_conductivity, eno_ns_surfref,
  enc_ncc_clcref, clutter_height, clutter_density, delta_h_diff,
  frq_mhz: double;
  const radio_climate, pol: integer;
  const pctTime, pctLoc, pctConf: double): double;

// other models
type
  TEnvironment = (evUrban, evSuburban, evRural);

var
  strEnvironment: array[TEnvironment] of
  string = ('City / Conservative', 'Suburban / Average', 'Rural / Optimistic');

// http://www.cl.cam.ac.uk/research/dtg/lce-pub/public/vsa23/VTC05_Empirical.pdf
function COST231pathLoss(f, TxH, RxH, d: double; const mode: TEnvironment): double;
function ECC33pathLoss(f, TxH, RxH, d: double; const mode: TEnvironment): double;
function EgliPathLoss(f, h1, h2, d: double): double;
function EricssonpathLoss(f, TxH, RxH, d: double; const Mode: TEnvironment): double;
function FSPLpathLoss(f, d: double): double;
function HATApathLoss(f, h_B, h_M, d: double; mode: TEnvironment): double;
function PlaneEarthLoss(d, TxH, RxH: double): double;
function SoilPathLoss(f, d, terdic: double): double;
function SUIpathLoss(f, TxH, RxH, d: double; mode: TEnvironment): double;

implementation

uses Math, ucomplex;

const
  THIRD = (1.0 / 3.0);

type
  prop_type = record
    aref: double;
    dist: double;
    hg: array[0..1] of double;
    rch: array[0..1] of double;
    wn: double;
    dh: double;
    dhd: double;
    ens: double;
    encc: double;
    cch: double;
    cd: double;
    gme: double;
    zgndreal: double;
    zgndimag: double;
    he: array[0..1] of double;
    dl: array[0..1] of double;
    the: array[0..1] of double;
    tiw: double;
    ght: double;
    ghr: double;
    rph: double;
    hht: double;
    hhr: double;
    tgh: double;
    tsgh: double;
    thera: double;
    thenr: double;
    rpl: integer;
    kwx: integer;
    mdp: integer;
    ptx: integer;
    los: integer;
  end;

  propv_type = record
    sgc: double;
    lvar: integer;
    mdvar: integer;
    klim: integer;
  end;

  propa_type = record
    dlsa: double;
    dx: double;
    ael: double;
    ak1: double;
    ak2: double;
    aed: double;
    emd: double;
    aes: double;
    ems: double;
    dls: array[0..1] of double;
    dla: double;
    tha: double;
  end;

  static_type = record
    wd1, xd1, afo, qk, aht, xht: double;                          // adiff
    toh, toho, roh, roho, dto, dto1,
    dtro, dro, dro2, drto, dtr, dhh1, dhh2, dtof, dto1f,
    drof, dro2f: double;                                          // addif2

    ad, rr, etq, h0s: double;                                     // ascat

    wls: double;                                                  // alos

    wlos, wscat: boolean;                                         // lrprop 2
    dmin, xae: double;                                            // lrprop 2

    kdv: integer;                                                 // avar
    dexa, de, vmd, vs0, sgl, sgtm, sgtp, sgtd, tgtd,              // avar
    gm, gp, cv1, cv2, yv1, yv2, yv3, csm1, csm2, ysm1, ysm2,
    ysm3, csp1, csp2, ysp1, ysp2, ysp3, csd1, zd, cfm1, cfm2,
    cfm3, cfp1, cfp2, cfp3: double;
    ws, w1: boolean;                                              // avar
  end;

function mymin(const i, j: integer): integer; inline; overload;
begin
  if (i < j) then
    Result := i
  else
    Result := j;
end;

function mymax(const i, j: integer): integer; inline; overload;
begin
  if (i > j) then
    Result := i
  else
    Result := j;
end;

function mymin(const a, b: double): double; inline; overload;
begin
  if (a < b) then
    Result := a
  else
    Result := b;
end;

function mymax(const a, b: double): double; inline; overload;
begin
  if (a > b) then
    Result := a
  else
    Result := b;
end;

function FORTRAN_DIM(const x, y: double): double; inline;
begin
      (* This performs the FORTRAN DIM function.  Result is x-y
         if x is greater than y; otherwise result is 0.0 *)
  if (x > y) then
    Result := x - y
  else
    Result := 0.0;
end;

function complexi(const re, im: double): complex; inline;
begin
  Result.re := re;
  Result.im := im;
end;

function absc(const c: complex): double;
begin
  Result := sqrt(c.re * c.re + c.im * c.im);
end;

function aknfe(const v2: double): double;
begin
  if (v2 < 5.76) then
    Result := 6.02 + 9.11 * sqrt(v2) - 1.27 * v2
  else
    Result := 12.953 + 10 * log10(v2);
end;

function fht(const x, pk: double): double;
var
  w, fhtv: double;
begin
  if (x < 200.0) then
  begin
    w := -ln(pk);

    if (pk < 1.0e-5) or (x * w * w * w > 5495.0) then
    begin
      fhtv := -117.0;

      if (x > 1.0) then
        fhtv := 40.0 * log10(x) + fhtv;
    end
    else
      fhtv := 2.5e-5 * x * x / pk - 8.686 * w - 15.0;
  end
  else
  begin
    fhtv := 0.05751 * x - 10.0 * log10(x);

    if (x < 2000.0) then
    begin
      w := 0.0134 * x * exp(-0.005 * x);
      fhtv := (1.0 - w) * fhtv + w * (40.0 * log10(x) - 117.0);
    end;
  end;
  Result := fhtv;
end;

function h0f(const r, et: double): double;
const
  a: array[0..4] of double = (25.0, 80.0, 177.0, 395.0, 705.0);
  b: array[0..4] of double = (24.0, 45.0, 68.0, 80.0, 105.0);
var
  q, x: double;
  h0fv, temp: double;
  it: integer;
begin
  it := trunc(et);

  if (it <= 0) then
  begin
    it := 1;
    q := 0.0;
  end
  else if (it >= 5) then
  begin
    it := 5;
    q := 0.0;
  end
  else
    q := et - it;

  (* x:=power(1.0/r,2.0); *)

  temp := 1.0 / r;
  x := temp * temp;

  h0fv := 4.343 * ln((a[it - 1] * x + b[it - 1]) * x + 1.0);

  if (q <> 0.0) then
    h0fv := (1.0 - q) * h0fv + q * 4.343 * ln((a[it] * x + b[it]) * x + 1.0);

  Result := h0fv;
end;

function ahd(td: double): double;
const
  a: array[0..2] of double = (133.4, 104.6, 71.8);
  b: array[0..2] of double = (0.332e-3, 0.212e-3, 0.157e-3);
  c: array[0..2] of double = (-4.343, -1.086, 2.171);
var
  i: integer;
begin
  if (td <= 10e3) then
    i := 0
  else if (td <= 70e3) then
    i := 1
  else
    i := 2;

  Result := a[i] + b[i] * td + c[i] * ln(td);
end;

function abq_alos(const r: complex): double;
begin
  Result := r.re * r.re + r.im * r.im;
end;

function saalos(d: double; const prop: prop_type; const {%H-}propa: propa_type): double;
var
  ensa, encca, q, dp, dx, tde, hc, ucrpc, ctip, tip, tic, stic, ctic, sta: double;
  ttc, cttc, crpc, ssnps, d1a, rsp, tsp, arte, zi, pd, pdk, hone, tvsr: double;
  saalosv: double;
  j: integer;
begin
  saalosv := 0.0;
  q := 0.0;
  if (d = 0.0) then
  begin
    tsp := 1.0;
    rsp := 0.0;
    d1a := 50.0;
    saalosv := 0.0;
  end
  else if (prop.hg[1] > prop.cch) then
  begin
    saalosv := 0.0;
  end
  else
  begin
    pd := d;
    pdk := pd / 1000.0;
    tsp := 1.0;
    rsp := 0.0;
    d1a := pd;
    (* at first, hone is transmitter antenna height
       relative to receive site ground level. *)
    hone := prop.tgh + prop.tsgh - (prop.rch[1] - prop.hg[1]);

    if (prop.tgh > prop.cch) then  (* for TX ant above all clutter height *)
    begin
      ensa := 1 + prop.ens * 0.000001;
      encca := 1 + prop.encc * 0.000001;
      dp := pd;

      for j := 0 to 4 do
      begin
        tde := dp / 6378137.0;
        hc := (prop.cch + 6378137.0) * (1 - cos(tde));
        dx := (prop.cch + 6378137.0) * sin(tde);
        ucrpc :=
          sqrt((hone - prop.cch + hc) * (hone - prop.cch + hc) +
          (dx * dx));
        ctip := (hone - prop.cch + hc) / ucrpc;
        tip := arccos(ctip);
        tic := tip + tde;
        tic := mymax(0.0, tic);
        stic := sin(tic);
        sta := (ensa / encca) * stic;
        ttc := arcsin(sta);
        cttc := sqrt(1 - (sin(ttc)) * (sin(ttc)));
        crpc := (prop.cch - prop.hg[1]) / cttc;
        if (crpc >= dp) then
        begin
          crpc := dp - 1 / dp;
        end;

        ssnps := (3.1415926535897 / 2) - tic;
        d1a := (crpc * sin(ttc)) / (1 - 1 / 6378137.0);
        dp := pd - d1a;
      end;

      ctic := cos(tic);

      (* if the ucrpc path touches the canopy before reaching the
         end of the ucrpc, the entry point moves toward the
         transmitter, extending the crpc and d1a. Estimating the d1a: *)

      if (ssnps <= 0.0) then
      begin
        d1a := mymin(0.1 * pd, 600.0);
        crpc := d1a;
        (* hone must be redefined as being barely above
           the canopy height with respect to the receiver
           canopy height, which despite the earth curvature
           is at or above the transmitter antenna height. *)
        hone := prop.cch + 1;
        rsp := 0.997;
        tsp := 1 - rsp;
      end
      else
      begin

        if (prop.ptx >= 1) then  (* polarity ptx is vertical or circular *)
        begin
          q := ((ensa * cttc - encca * ctic) /
            (ensa * cttc + encca * ctic));
          rsp := q * q;
          tsp := 1 - rsp;

          if (prop.ptx = 2) then   (* polarity is circular - new *)
          begin
            q := ((ensa * ctic - encca * cttc) /
              (ensa * ctic + encca * cttc));
            rsp :=
              ((ensa * cttc - encca * ctic) /
              (ensa * cttc + encca * ctic));
            rsp := (q * q + rsp * rsp) / 2;
            tsp := 1 - rsp;
          end;
        end
        else   (* ptx is 0, horizontal, or undefined *)
        begin
          q := ((ensa * ctic - encca * cttc) /
            (ensa * ctic + encca * cttc));
          rsp := q * q;
          tsp := 1 - rsp;
        end;
      end;
      (* tvsr is defined as tx ant height above receiver ant height *)
      tvsr := mymax(0.0, prop.tgh + prop.tsgh - prop.rch[1]);

      if (d1a < 50.0) then
      begin
        arte := 0.0195 * crpc - 20 * log10(tsp);
      end
      else
      begin
        if (d1a < 225.0) then
        begin

          if (tvsr > 1000.0) then
          begin
            q := d1a * (0.03 * exp(-0.14 * pdk));
          end
          else
          begin
            q := d1a * (0.07 * exp(-0.17 * pdk));
          end;

          arte :=
            q + (0.7 * pdk - mymax(0.01, log10(prop.wn * 47.7) -
            2)) * (prop.hg[1] / hone);
        end
        else
        begin
          q := 0.00055 * (pdk) + log10(pdk) *
            (0.041 - 0.0017 * sqrt(hone) + 0.019);

          arte :=
            d1a * q - (18 * log10(rsp)) / (exp(hone / 37.5));

          zi := 1.5 * sqrt(hone - prop.cch);

          if (pdk > zi) then
          begin
            q := (pdk - zi) * 10.2 *
              ((sqrt(mymax(0.01, log10(prop.wn * 47.7) -
              2.0))) / (100 - zi));
          end
          else
          begin
            q := ((zi - pdk) / zi) *
              (-20.0 * mymax(0.01, log10(prop.wn *
              47.7) - 2.0)) / sqrt(hone);
          end;
          arte := arte + q;

        end;
      end;
    end
    else   (* for TX at or below clutter height *)
    begin
      q := (prop.cch - prop.tgh) * (2.06943 - 1.56184 *
        exp(1 / prop.cch - prop.tgh));
      q := q + (17.98 - 0.84224 * (prop.cch - prop.tgh)) *
        exp(-0.00000061 * pd);
      arte := q + 1.34795 * 20 * log10(pd + 1.0);
      arte :=
        arte - (mymax(0.01, log10(prop.wn * 47.7) - 2)) *
        (prop.hg[1] / prop.tgh);
    end;
    saalosv := arte;
  end;
  Result := saalosv;
end;

function adiff(const d: double; const prop: prop_type; const propa: propa_type;
  var static: static_type): double;
var
  prop_zgnd: complex;
  a, q, pk, ds, th, wa, ar, wd, adiffv: double;
  j: integer;
begin
  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);

  if (d = 0) then
  begin
    q := prop.hg[0] * prop.hg[1];
    static.qk := prop.he[0] * prop.he[1] - q;

    if (prop.mdp < 0.0) then
      q := q + 10.0;

    static.wd1 := sqrt(1.0 + static.qk / q);
    static.xd1 := propa.dla + propa.tha / prop.gme;
    q := (1.0 - 0.8 * exp(-propa.dlsa / 50e3)) * prop.dh;
    q := q * (0.78 * exp(-power(q / 16.0, 0.25)));
    static.afo :=
      mymin(15.0, 2.171 * ln(1.0 + 4.77e-4 * prop.hg[0] *
      prop.hg[1] * prop.wn * q));
    static.qk := 1.0 / absc(prop_zgnd);
    static.aht := 20.0;
    static.xht := 0.0;

    for j := 0 to 1 do
    begin
      (* a=0.5*power(prop.dl[j],2.0)/prop.he[j]; *)
      a := 0.5 * (prop.dl[j] * prop.dl[j]) / prop.he[j];
      wa := power(a * prop.wn, THIRD);
      pk := static.qk / wa;
      q := (1.607 - pk) * 151.0 * wa * prop.dl[j] / a;
      static.xht := static.xht + q;
      static.aht := static.aht + fht(q, pk);
    end;

    adiffv := 0.0;
  end
  else
  begin
    th := propa.tha + d * prop.gme;
    ds := d - propa.dla;
    (* q=0.0795775*prop.wn*ds*power(th,2.0); *)
    q := 0.0795775 * prop.wn * ds * th * th;
    adiffv :=
      aknfe(q * prop.dl[0] / (ds + prop.dl[0])) +
      aknfe(q * prop.dl[1] / (ds + prop.dl[1]));
    a := ds / th;
    wa := power(a * prop.wn, THIRD);
    pk := static.qk / wa;
    q := (1.607 - pk) * 151.0 * wa * th + static.xht;
    ar := 0.05751 * q - 4.343 * ln(q) - static.aht;
    q := (static.wd1 + static.xd1 / d) * mymin(
      ((1.0 - 0.8 * exp(-d / 50e3)) * prop.dh * prop.wn), 6283.2);
    wd := 25.1 / (25.1 + sqrt(q));
    adiffv := ar * wd + (1.0 - wd) * adiffv + static.afo;
  end;

  Result := adiffv;
end;

function adiff2(const d: double; var prop: prop_type; const propa: propa_type;
  var static: static_type): double;
var
  prop_zgnd: complex;
  a, q, pk, rd, ds, dsl, th, wa: double;
  sf2, vv, kedr, arp, sdr, pd, srp, kem, csd, sdl, adiffv2, closs: double;
begin
  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);
  kedr := 0.0;
  arp := 0.0;
  sdr := 0.0;
  pd := 0.0;
  srp := 0.0;
  kem := 0.0;
  csd := 0.0;
  sdl := 0.0;
  adiffv2 := 0.0;
  closs := 0.0;

  (* sf1=1.0; *)(* average empirical hilltop foliage scatter factor for 1 obstruction  *)
  sf2 := 1.0;    (* average empirical hilltop foliage scatter factor for 2 obstructions *)

  (* dfdh=prop.dh; *)
  (* ec=0.5*prop.gme; *)

  (* adiff2 must first be run with d==0.0 to set up coefficients *)
  if (d = 0) then
  begin
    q := prop.hg[0] * prop.hg[1];
    static.qk := prop.he[0] * prop.he[1] - q;
    (* dhec=2.73; *)

    if (prop.mdp < 0.0) then
      q := q + 10.0;

    (* coefficients for a standard four radii, rounded earth computation are prepared *)
    static.wd1 := sqrt(1.0 + static.qk / q);
    static.xd1 := propa.dla + propa.tha / prop.gme;
    q := (1.0 - 0.8 * exp(-propa.dlsa / 50e3)) * prop.dh;
    q := q * (0.78 * exp(-power(q / 16.0, 0.25)));
    static.qk := 1.0 / absc(prop_zgnd);
    static.aht := 20.0;
    static.xht := 0.0;
    a := 0.5 * (prop.dl[0] * prop.dl[0]) / prop.he[0];
    wa := power(a * prop.wn, THIRD);
    pk := static.qk / wa;
    q := (1.607 - pk) * 151.0 * wa * prop.dl[0] / a;
    static.xht := q;
    static.aht := static.aht + fht(q, pk);

    if ((trunc(prop.dl[1]) = 0.0) or (prop.the[1] > 0.2)) then
    begin
      static.xht := static.xht + static.xht;
      static.aht := static.aht + (static.aht - 20.0);
    end
    else
    begin
      a := 0.5 * (prop.dl[1] * prop.dl[1]) / prop.he[1];
      wa := power(a * prop.wn, THIRD);
      pk := static.qk / wa;
      q := (1.607 - pk) * 151.0 * wa * prop.dl[1] / a;
      static.xht := static.xht + q;
      static.aht := static.aht + fht(q, pk);
    end;
    adiffv2 := 0.0;
  end
  else
  begin
    th := propa.tha + d * prop.gme;

    dsl := mymax(d - propa.dla, 0.0);
    ds := d - propa.dla;
    a := ds / th;
    wa := power(a * prop.wn, THIRD);
    pk := static.qk / wa;
    static.toh :=
      prop.hht - (prop.rch[0] - prop.dl[0] *
      ((prop.rch[1] - prop.rch[0]) / prop.dist));
    static.roh :=
      prop.hhr - (prop.rch[0] - (prop.dist - prop.dl[1]) *
      ((prop.rch[1] - prop.rch[0]) / prop.dist));
    static.toho :=
      prop.hht - (prop.rch[0] - (prop.dl[0] + dsl) *
      ((prop.hhr - prop.rch[0]) / (prop.dist - prop.dl[1])));
    static.roho :=
      prop.hhr - (prop.hht - dsl * ((prop.rch[1] - prop.hht) / dsl));
    static.dto := sqrt(prop.dl[0] * prop.dl[0] + static.toh * static.toh);
    static.dto := static.dto + (prop.gme * prop.dl[0]);
    static.dto1 := sqrt(prop.dl[0] * prop.dl[0] + static.toho * static.toho);
    static.dto1 := static.dto1 + (prop.gme * prop.dl[0]);
    static.dtro :=
      sqrt((prop.dl[0] + dsl) * (prop.dl[0] + dsl) + prop.hhr * prop.hhr);
    static.dtro := static.dtro + (prop.gme * (prop.dl[0] + dsl));
    static.drto :=
      sqrt((prop.dl[1] + dsl) * (prop.dl[1] + dsl) + prop.hht * prop.hht);
    static.drto := static.drto + (prop.gme * (prop.dl[1] + dsl));
    static.dro := sqrt(prop.dl[1] * prop.dl[1] + static.roh * static.roh);
    static.dro := static.dro + (prop.gme * (prop.dl[1]));
    static.dro2 := sqrt(prop.dl[1] * prop.dl[1] + static.roho * static.roho);
    static.dro2 := static.dro2 + (prop.gme * (prop.dl[1]));
    static.dtr :=
      sqrt(prop.dist * prop.dist + (prop.rch[0] - prop.rch[1]) *
      (prop.rch[0] - prop.rch[1]));
    static.dtr := static.dtr + (prop.gme * prop.dist);
    static.dhh1 :=
      sqrt((prop.dist - propa.dla) * (prop.dist - propa.dla) +
      static.toho * static.toho);
    static.dhh1 := static.dhh1 + (prop.gme * (prop.dist - propa.dla));
    static.dhh2 :=
      sqrt((prop.dist - propa.dla) * (prop.dist - propa.dla) +
      static.roho * static.roho);
    static.dhh2 := static.dhh2 + (prop.gme * (prop.dist - propa.dla));

    (* for 1 obst tree base path *)
    static.dtof :=
      sqrt(prop.dl[0] * prop.dl[0] + (static.toh - prop.cch) *
      (static.toh - prop.cch));
    static.dtof := static.dtof + (prop.gme * prop.dl[0]);
    static.dto1f :=
      sqrt(prop.dl[0] * prop.dl[0] + (static.toho - prop.cch) *
      (static.toho - prop.cch));
    static.dto1f := static.dto1f + (prop.gme * prop.dl[0]);
    static.drof :=
      sqrt(prop.dl[1] * prop.dl[1] + (static.roh - prop.cch) *
      (static.roh - prop.cch));
    static.drof := static.drof + (prop.gme * (prop.dl[1]));
    static.dro2f :=
      sqrt(prop.dl[1] * prop.dl[1] + (static.roho - prop.cch) *
      (static.roho - prop.cch));
    static.dro2f := static.dro2f + (prop.gme * (prop.dl[1]));

    (* saalos coefficients preset for post-obstacle receive path *)
    prop.tgh := prop.cch + 1.0;
    prop.tsgh := prop.hhr;
    rd := prop.dl[1];

    (* two obstacle diffraction calculation *)
    if (trunc(ds) > 0) then  (* there are 2 obstacles *)
    begin
      if (trunc(prop.dl[1]) > 0.0) then  (* receive site past 2nd peak *)
      begin
        (* rounding attenuation *)
        q := (1.607 - pk) * 151.0 * wa * th + static.xht;
        (* ar=0.05751*q-10*log10(q)-aht; *)

        (* knife edge vs round weighting *)
        q := (1.0 - 0.8 * exp(-d / 50e3)) * prop.dh;
        q := (static.wd1 + static.xd1 / d) * mymin((q * prop.wn), 6283.2);
        (* wd=25.1/(25.1+sqrt(q)); *)

        q := 0.6365 * prop.wn;

        if (prop.the[1] < 0.2) then  (* receive grazing angle below 0.2 rad *)
        begin
          (* knife edge attenuation for two obstructions *)

          if (prop.hht < 3400) then   (* if below tree line, foliage top loss *)
          begin
            vv := q * abs(static.dto1 + static.dhh1 - static.dtro);
            adiffv2 := -18.0 + sf2 * aknfe(vv);
          end
          else
          begin
            vv := q * abs(static.dto1 + static.dhh1 - static.dtro);
            adiffv2 := aknfe(vv);
          end;

          if (prop.hhr < 3400) then
          begin
            vv := q * abs(static.dro2 + static.dhh2 - static.drto);
            adiffv2 := adiffv2 + (-18.0 + sf2 * aknfe(vv));
          end
          else
          begin
            vv := q * abs(static.dro2 + static.dhh2 - static.drto);
            adiffv2 := adiffv2 + aknfe(vv);
          end;
          (* finally, add clutter loss *)
          closs := saalos(rd, prop, propa);
          adiffv2 := adiffv2 + mymin(22.0, closs);

        end
        else (* rcvr site too close to 2nd obs *)
        begin
          (* knife edge attenuation for 1st obs *)

          if (prop.hht < 3400) then
          begin
            vv := q * abs(static.dto1 + static.dhh1 - static.dtro);
            adiffv2 := -18.0 + sf2 * aknfe(vv);
          end
          else
          begin
            vv := q * abs(static.dto1 + static.dhh1 - static.dtro);
            adiffv2 := aknfe(vv);
          end;

          (* weighted calc. of knife vs rounded edge
             adiffv2=ar*wd+(1.0-wd)*adiffv2; *)

          (* clutter path loss past 2nd peak *)
          if (prop.the[1] < 1.22) then
          begin
            rd := prop.dl[1];

            if (prop.the[1] > 0.6) then  (* through foliage downhill *)
            begin
              prop.tgh := prop.cch;
            end
            else
              (* close to foliage, rcvr in foliage downslope *)
            begin

              vv := 0.6365 * prop.wn * abs(static.dro2 + static.dhh2 -
                static.drto);
            end;
            adiffv2 := adiffv2 + aknfe(vv);
            closs := saalos(rd, prop, propa);
            adiffv2 := adiffv2 + mymin(closs, 22.0);
          end
          else  (* rcvr very close to bare cliff or skyscraper *)
          begin

            adiffv2 := 5.8 + 25.0;
          end;
        end;
      end
      else (* receive site is atop a 2nd peak *)
      begin
        vv := 0.6365 * prop.wn * abs(static.dto + static.dro - static.dtr);
        adiffv2 := 5.8 + aknfe(vv);
      end;
    end
    else (* for single obstacle *)
    begin

      if (trunc(prop.dl[1]) > 0.0) then  (* receive site past 1st peak *)
      begin

        if (prop.the[1] < 0.2) then  (* receive grazing angle less than .2 radians *)
        begin
          vv := 0.6365 * prop.wn * abs(static.dto + static.dro - static.dtr);

          if (prop.hht < 3400) then
          begin
            sdl := 18.0;
            sdl := power(10, (-sdl / 20));
            (* ke phase difference with respect to direct t-r line *)
            kedr :=
              0.159155 * prop.wn * abs(static.dto + static.dro - static.dtr);
            arp := abs(kedr - (int(kedr)));
            kem := aknfe(vv);
            kem := power(10, (-kem / 20));
            (* scatter path phase with respect to direct t-r line *)
            sdr :=
              0.5 + 0.159155 * prop.wn * abs(static.dtof +
              static.drof - static.dtr);
            srp := abs(sdr - (trunc(sdr)));
            (* difference between scatter and ke phase in radians *)
            pd := 6.283185307 * abs(srp - arp);
            (* report pd prior to restriction
               keep pd between 0 and pi radians and adjust for 3&4 quadrant *)
            if (pd >= 3.141592654) then
            begin
              pd := 6.283185307 - pd;
              csd :=
                abq_alos(complexi(sdl, 0) + complexi(kem *
                -cos(pd), kem * -sin(pd)));
            end
            else
            begin
              csd :=
                abq_alos(complexi(sdl, 0) + complexi(kem *
                cos(pd), kem * sin(pd)));
            end;
            (*csd=mymax(csd,0.0009); limits maximum loss value to 30.45 db *)
            adiffv2 := -3.71 - 10 * log10(csd);
          end
          else
          begin
            adiffv2 := aknfe(vv);
          end;
          (* finally, add clutter loss *)
          closs := saalos(rd, prop, propa);
          adiffv2 := adiffv2 + mymin(closs, 22.0);
        end
        else  (* receive grazing angle too high *)
        begin

          if (prop.the[1] < 1.22) then
          begin
            rd := prop.dl[1];

            if (prop.the[1] > 0.6) then  (* through foliage downhill *)
            begin
              prop.tgh := prop.cch;
            end
            else   (* downhill slope just above foliage  *)
            begin

              vv := 0.6365 * prop.wn * abs(static.dto + static.dro -
                static.dtr);
              adiffv2 := aknfe(vv);
            end;
            closs := saalos(rd, prop, propa);
            adiffv2 := adiffv2 + mymin(22.0, closs);
          end
          else   (* receiver very close to bare cliff or skyscraper *)
          begin
            adiffv2 := 5.8 + 25.0;
          end;
        end;
      end
      else   (* if occurs, receive site atop first peak  *)
      begin
        adiffv2 := 5.8;
      end;
    end;
  end;
  Result := adiffv2;
end;

function ascat(const d: double; const prop: prop_type; propa: propa_type;
  var static: static_type): double;
var
  h0, r1, r2, z0, ss, et, ett, th, q: double;
  ascatv, temp: double;
begin
  if (d = 0.0) then
  begin
    static.ad := prop.dl[0] - prop.dl[1];
    static.rr := prop.he[1] / prop.rch[0];

    if (static.ad < 0.0) then
    begin
      static.ad := -static.ad;
      static.rr := 1.0 / static.rr;
    end;

    static.etq := (5.67e-6 * prop.ens - 2.32e-3) * prop.ens + 0.031;
    static.h0s := -15.0;
    ascatv := 0.0;
  end
  else
  begin
    if (static.h0s > 15.0) then
    begin
      h0 := static.h0s;
    end
    else
    begin
      th := prop.the[0] + prop.the[1] + d * prop.gme;
      r2 := 2.0 * prop.wn * th;
      r1 := r2 * prop.he[0];
      r2 := r2 * prop.he[1];

      if (r1 < 0.2) and (r2 < 0.2) then
      begin
        Result := 1001.0;  // <==== early return
        exit;
      end;

      ss := (d - static.ad) / (d + static.ad);
      q := static.rr / ss;
      ss := mymax(0.1, ss);
      q := mymin(mymax(0.1, q), 10.0);
      z0 := (d - static.ad) * (d + static.ad) * th * 0.25 / d;
      (* et:=(etq*exp(-power(mymin(1.7,z0/8.0e3),6.0))+1.0)*z0/1.7556e3; *)

      temp := mymin(1.7, z0 / 8.0e3);
      temp := temp * temp * temp * temp * temp * temp;
      et := (static.etq * exp(-temp) + 1.0) * z0 / 1.7556e3;

      ett := mymax(et, 1.0);
      h0 := (h0f(r1, ett) + h0f(r2, ett)) * 0.5;
      h0 := h0 + mymin(h0, (1.38 - ln(ett)) * ln(ss) * ln(q) * 0.49);
      h0 := FORTRAN_DIM(h0, 0.0);

      if (et < 1.0) then
      begin
        (* h0=et*h0+(1.0-et)*4.343*ln(power((1.0+1.4142/r1)*(1.0+1.4142/r2),2.0)*(r1+r2)/(r1+r2+2.8284)); *)

        temp :=
          ((1.0 + 1.4142 / r1) * (1.0 + 1.4142 / r2));
        h0 := et * h0 + (1.0 - et) * 4.343 * ln(
          (temp * temp) * (r1 + r2) /
          (r1 + r2 + 2.8284));
      end;

      if (h0 > 15.0) and (static.h0s >= 0.0) then
      begin
        h0 := static.h0s;
      end;
    end;

    static.h0s := h0;
    th := propa.tha + d * prop.gme;
    (* ascatv=ahd(th*d)+4.343*ln(47.7*prop.wn*power(th,4.0))-0.1*(prop.ens-301.0)*exp(-th*d/40e3)+h0; *)
    ascatv :=
      ahd(th * d) + 4.343 * ln(47.7 * prop.wn * (th * th * th * th)) -
      0.1 * (prop.ens - 301.0) * exp(-th * d / 40e3) + h0;
  end;

  Result := ascatv;
end;

function qerfi(const q: double): double;
const
  c0 = 2.515516698;
  c1 = 0.802853;
  c2 = 0.010328;
  d1 = 1.432788;
  d2 = 0.189269;
  d3 = 0.001308;

var
  x, t, v: double;
begin
  x := 0.5 - q;
  t := mymax(0.5 - abs(x), 0.000001);
  t := sqrt(-2.0 * ln(t));
  v := t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1.0);

  if (x < 0.0) then
    v := -v;

  Result := v;
end;

procedure qlrps(const fmhz, zsys, en0: double; const ipol: integer;
  const eps, sgm: double; var prop: prop_type);
var
  gma: double;
  zq, prop_zgnd: complex;
begin
  gma := 157e-9;

  prop.wn := fmhz / 47.7;
  prop.ens := en0;

  if (zsys <> 0.0) then
    prop.ens := prop.ens * exp(-zsys / 9460.0);

  prop.gme := gma * (1.0 - 0.04665 * exp(prop.ens / 179.3));
  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);
  zq := complexi(eps, 376.62 * sgm / prop.wn);
  prop_zgnd := csqrt(zq - 1.0);

  if (ipol <> 0.0) then
    prop_zgnd := prop_zgnd / zq;

  prop.zgndreal := prop_zgnd.re;
  prop.zgndimag := prop_zgnd.im;
end;

function alos(const d: double; const prop: prop_type;
  const propa: propa_type; var static: static_type): double;
var
  r, prop_zgnd: complex;
  s, sps, q: double;
  alosv: double;
begin
  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);

  if (d = 0.0) then
  begin
    static.wls :=
      0.021 / (0.021 + prop.wn * prop.dh / mymax(10e3, propa.dlsa));
    alosv := 0.0;
  end
  else
  begin
    q := (1.0 - 0.8 * exp(-d / 50e3)) * prop.dh;
    s := 0.78 * q * exp(-power(q / 16.0, 0.25));
    q := prop.he[0] + prop.he[1];
    sps := q / sqrt(d * d + q * q);
    r := (sps - prop_zgnd) / (sps + prop_zgnd) *
      exp(-mymin(10.0, prop.wn * s * sps));
    q := abq_alos(r);

    if (q < 0.25) or (q < sps) then
      r := r * sqrt(sps / q);

    alosv := propa.emd * d + propa.aed;
    q := prop.wn * prop.he[0] * prop.he[1] * 2.0 / d;

    if (q > 1.57) then
      q := 3.14 - 2.4649 / q;

    alosv :=
      (-4.343 * ln(abq_alos(complexi(cos(q), -sin(q)) + r)) - alosv) *
      static.wls + alosv;
  end;
  Result := alosv;
end;

function alos2(const d: double; var prop: prop_type; const propa: propa_type): double;
var
  r, prop_zgnd: complex;
  cd, cr, dr, hr, hrg, ht, htg, hrp, re, s, sps, q, pd, drh: double;
  alosv: double;
begin
  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);

  cd := 0.0;
  cr := 0.0;
  htg := prop.hg[0];
  hrg := prop.hg[1];
  ht := prop.ght;
  hr := prop.ghr;
  (* rp:=prop.rpl; *)
  hrp := prop.rph;
  pd := prop.dist;

  if (d = 0.0) then
  begin
    alosv := 0.0;
  end
  else
  begin
    q := prop.he[0] + prop.he[1];
    sps := q / sqrt(pd * pd + q * q);
    q := (1.0 - 0.8 * exp(-pd / 50e3)) * prop.dh;

    if (prop.mdp < 0) then
    begin
      dr := pd / (1 + hrg / htg);

      if (dr < (0.5 * pd)) then
      begin
        drh :=
          6378137.0 - sqrt(-(0.5 * pd) * (0.5 * pd) + 6378137.0 *
          6378137.0 + (0.5 * pd - dr) * (0.5 * pd - dr));
      end
      else
      begin
        drh :=
          6378137.0 - sqrt(-(0.5 * pd) * (0.5 * pd) + 6378137.0 *
          6378137.0 + (dr - 0.5 * pd) * (dr - 0.5 * pd));
      end;

      if ((sps < 0.05) and (prop.cch > hrg) and (prop.dist < prop.dl[0])) then
        (* if far from transmitter and receiver below canopy *)
      begin
        cd :=
          mymax(0.01, pd * (prop.cch - hrg) / (htg - hrg));
        cr := mymax(0.01, pd - dr + dr * (prop.cch - drh) / htg);
        q := ((1.0 - 0.8 * exp(-pd / 50e3)) * prop.dh *
          (mymin(-20 * log10(cd / cr), 1.0)));
      end;
    end;

    s := 0.78 * q * exp(-power(q / 16.0, 0.25));
    q := exp(-mymin(10.0, prop.wn * s * sps));
    r := q * (sps - prop_zgnd) / (sps + prop_zgnd);
    q := abq_alos(r);
    q := mymin(q, 1.0);

    if (q < 0.25) or (q < sps) then
    begin
      r := r * sqrt(sps / q);
    end;
    q := prop.wn * prop.he[0] * prop.he[1] / (pd * 3.1415926535897);

    if (prop.mdp < 0) then
    begin
      q := prop.wn * ((ht - hrp) * (hr - hrp)) / (pd * 3.1415926535897);
    end;
    q := q - floor(q);

    if (q < 0.5) then
    begin
      q := q * 3.1415926535897;
    end
    else
    begin
      q := (1 - q) * 3.1415926535897;
    end;

    (* no longer valid complex conjugate removed
       by removing minus sign from in front of sin function *)
    re := abq_alos(complexi(cos(q), sin(q)) + r);
    alosv := -10 * log10(re);
    prop.tgh := prop.hg[0];  (*tx above gnd hgt set to antenna height AGL *)
    prop.tsgh := prop.rch[0] - prop.hg[0];  (* tsgh set to tx site gl AMSL *)

    if ((prop.hg[1] < prop.cch) and (prop.thera < 0.785) and
      (prop.thenr < 0.785)) then
    begin
      if (sps < 0.05) then
      begin
        alosv := alosv + saalos(pd, prop, propa);
      end
      else
      begin
        alosv := saalos(pd, prop, propa);
      end;
    end;
  end;
  alosv := mymin(22.0, alosv);
  Result := alosv;
end;

procedure qlra(const kst: array of integer; const klimx, mdvarx: integer;
  var prop: prop_type; var propv: propv_type);
var
  q: double;
  j: integer;
begin
  for j := 0 to 1 do
  begin
    if (kst[j] <= 0) then
    begin
      prop.he[j] := prop.hg[j];
    end
    else
    begin
      q := 4.0;

      if (kst[j] <> 1) then
        q := 9.0;

      if (prop.hg[j] < 5.0) then
        q := q * sin(0.3141593 * prop.hg[j]);

      prop.he[j] :=
        prop.hg[j] + (1.0 + q) *
        exp(-mymin(20.0, 2.0 * prop.hg[j] / mymax(1e-3, prop.dh)));
    end;

    q := sqrt(2.0 * prop.he[j] / prop.gme);
    prop.dl[j] :=
      q * exp(-0.07 * sqrt(prop.dh / mymax(prop.he[j], 5.0)));
    prop.the[j] :=
      (0.65 * prop.dh * (q / prop.dl[j] - 1.0) - 2.0 * prop.he[j]) / q;
  end;

  prop.mdp := 1;
  propv.lvar := mymax(propv.lvar, 3);

  if (mdvarx >= 0) then
  begin
    propv.mdvar := mdvarx;
    propv.lvar := mymax(propv.lvar, 4);
  end;

  if (klimx > 0) then
  begin
    propv.klim := klimx;
    propv.lvar := 5;
  end;
end;

procedure lrprop(const d: double; var prop: prop_type; var propa: propa_type;
  var static: static_type);
var
  prop_zgnd: complex;
  a0, a1, a2, a3, a4, a5, a6: double;
  d0, d1, d2, d3, d4, d5, d6: double;
  wq: boolean;
  q: double;
  j: integer;

begin
  (* PaulM_lrprop used for ITM *)

  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);

  if (prop.mdp <> 0) then
  begin
    for j := 0 to 1 do
      propa.dls[j] := sqrt(2.0 * prop.he[j] / prop.gme);

    propa.dlsa := propa.dls[0] + propa.dls[1];
    propa.dla := prop.dl[0] + prop.dl[1];
    propa.tha :=
      mymax(prop.the[0] + prop.the[1], -propa.dla * prop.gme);
    static.wlos := False;
    static.wscat := False;

    if (prop.wn < 0.838) or (prop.wn > 210.0) then
      prop.kwx := mymax(prop.kwx, 1);

    for j := 0 to 1 do
      if (prop.hg[j] < 1.0) or (prop.hg[j] > 1000.0) then
        prop.kwx := mymax(prop.kwx, 1);

    for j := 0 to 1 do
      if ((abs(prop.the[j]) > 200e-3) or (prop.dl[j] < 0.1 * propa.dls[j]) or
        (prop.dl[j] > 3.0 * propa.dls[j])) then
        prop.kwx := mymax(prop.kwx, 3);

    if (prop.ens < 250.0) or (prop.ens > 400.0) or (prop.gme < 75e-9) or
      (prop.gme > 250e-9) or (prop_zgnd.re <= abs(prop_zgnd.im)) or
      (prop.wn < 0.419) or (prop.wn > 420.0) then
      prop.kwx := 4;

    for j := 0 to 1 do
      if (prop.hg[j] < 0.5) or (prop.hg[j] > 3000.0) then
        prop.kwx := 4;

    static.dmin := abs(prop.he[0] - prop.he[1]) / 200e-3;
    q := adiff(0.0, prop, propa, static);
    (* xae := power(prop.wn*power(prop.gme,2.),-THIRD); -- JDM made argument 2 a double *)
    static.xae := power(prop.wn * (prop.gme * prop.gme), -THIRD);  (* No 2nd power() *)
    d3 := mymax(propa.dlsa, 1.3787 * static.xae + propa.dla);
    d4 := d3 + 2.7574 * static.xae;
    a3 := adiff(d3, prop, propa, static);
    a4 := adiff(d4, prop, propa, static);
    propa.emd := (a4 - a3) / (d4 - d3);
    propa.aed := a3 - propa.emd * d3;
  end;

  if (prop.mdp >= 0) then
  begin
    prop.mdp := 0;
    prop.dist := d;
  end;

  if (prop.dist > 0.0) then
  begin
    if (prop.dist > 1000e3) then
      prop.kwx := mymax(prop.kwx, 1);

    if (prop.dist < static.dmin) then
      prop.kwx := mymax(prop.kwx, 3);

    if (prop.dist < 1e3) or (prop.dist > 2000e3) then
      prop.kwx := 4;
  end;

  if (prop.dist < propa.dlsa) then
  begin
    if (not static.wlos) then
    begin
      q := alos(0.0, prop, propa, static);
      d2 := propa.dlsa;
      a2 := propa.aed + d2 * propa.emd;
      d0 := 1.908 * prop.wn * prop.he[0] * prop.he[1];

      if (propa.aed >= 0.0) then
      begin
        d0 := mymin(d0, 0.5 * propa.dla);
        d1 := d0 + 0.25 * (propa.dla - d0);
      end
      else
        d1 := mymax(-propa.aed / propa.emd, 0.25 * propa.dla);

      a1 := alos(d1, prop, propa, static);
      wq := False;

      if (d0 < d1) then
      begin
        a0 := alos(d0, prop, propa, static);
        q := ln(d2 / d0);
        propa.ak2 :=
          mymax(0.0, ((d2 - d0) * (a1 - a0) - (d1 - d0) * (a2 - a0)) /
          ((d2 - d0) * ln(d1 / d0) -
          (d1 - d0) * q));
        wq := (propa.aed >= 0.0) or (propa.ak2 > 0.0);

        if (wq) then
        begin
          propa.ak1 :=
            (a2 - a0 - propa.ak2 * q) / (d2 - d0);

          if (propa.ak1 < 0.0) then
          begin
            propa.ak1 := 0.0;
            propa.ak2 :=
              FORTRAN_DIM(a2, a0) / q;

            if (propa.ak2 = 0.0) then
              propa.ak1 := propa.emd;
          end;
        end
        else
        begin
          propa.ak2 := 0.0;
          propa.ak1 := (a2 - a1) / (d2 - d1);

          if (propa.ak1 <= 0.0) then
            propa.ak1 := propa.emd;
        end;
      end
      else
      begin
        propa.ak1 := (a2 - a1) / (d2 - d1);
        propa.ak2 := 0.0;

        if (propa.ak1 <= 0.0) then
          propa.ak1 := propa.emd;
      end;

      propa.ael := a2 - propa.ak1 * d2 - propa.ak2 * ln(d2);
      static.wlos := True;
    end;

    if (prop.dist > 0.0) then
      prop.aref :=
        propa.ael + propa.ak1 * prop.dist + propa.ak2 * ln(prop.dist);
  end;

  if (prop.dist <= 0.0) or (prop.dist >= propa.dlsa) then
  begin
    if (not static.wscat) then
    begin
      q := ascat(0.0, prop, propa, static);
      d5 := propa.dla + 200e3;
      d6 := d5 + 200e3;
      a6 := ascat(d6, prop, propa, static);
      a5 := ascat(d5, prop, propa, static);

      if (a5 < 1000.0) then
      begin
        propa.ems := (a6 - a5) / 200e3;
        propa.dx :=
          mymax(propa.dlsa, mymax(propa.dla + 0.3 * static.xae *
          ln(47.7 * prop.wn), (a5 - propa.aed - propa.ems * d5) /
          (propa.emd - propa.ems)));
        propa.aes :=
          (propa.emd - propa.ems) * propa.dx + propa.aed;
      end
      else
      begin
        propa.ems := propa.emd;
        propa.aes := propa.aed;
        propa.dx := 10.e6;
      end;

      static.wscat := True;
    end;

    if (prop.dist > propa.dx) then
      prop.aref := propa.aes + propa.ems * prop.dist
    else
      prop.aref := propa.aed + propa.emd * prop.dist;
  end;

  prop.aref := mymax(prop.aref, 0.0);
end;

procedure lrprop2(const d: double; var prop: prop_type; var propa: propa_type;
  var static: static_type);
var
  prop_zgnd: complex;
  pd1: double;
  a0, a1, a2, a3, a4, a5, a6, iw: double;
  d0, d1, d2, d3, d4, d5, d6: double;
  wq: boolean;
  q: double;
  j: integer;

begin
  (* ITWOM_lrprop2 *)
  prop_zgnd := complexi(prop.zgndreal, prop.zgndimag);

  iw := prop.tiw;
  pd1 := prop.dist;
  propa.dx := 2000000.0;

  if (prop.mdp <> 0) then
  begin (* if oper. mode is not 0, i.e. not area mode ongoing *)
    for j := 0 to 1 do
      propa.dls[j] := sqrt(2.0 * prop.he[j] / prop.gme);

    propa.dlsa := propa.dls[0] + propa.dls[1];
    propa.dlsa := mymin(propa.dlsa, 1000000.0);
    propa.dla := prop.dl[0] + prop.dl[1];
    propa.tha :=
      mymax(prop.the[0] + prop.the[1], -propa.dla * prop.gme);
    static.wlos := False;
    static.wscat := False;

    (*checking for parameters-in-range, error codes set if not *)

    if (prop.wn < 0.838) or (prop.wn > 210.0) then
      prop.kwx := mymax(prop.kwx, 1);

    for j := 0 to 1 do
      if (prop.hg[j] < 1.0) or (prop.hg[j] > 1000.0) then
        prop.kwx := mymax(prop.kwx, 1);

    if (abs(prop.the[0]) > 200e-3) then
      prop.kwx := mymax(prop.kwx, 3);

    if (abs(prop.the[1]) > 1.220) then
      prop.kwx := mymax(prop.kwx, 3);

    (*for j := 0 to 1 do
       if (prop.dl[j]<0.1*propa.dls[j]) or  (prop.dl[j]>3.0*propa.dls[j]) then
        prop.kwx := mymax(prop.kwx,3);
                *)

    if (prop.ens < 250.0) or (prop.ens > 400.0) or (prop.gme < 75e-9) or
      (prop.gme > 250e-9) or (prop_zgnd.re <= abs(prop_zgnd.im)) or
      (prop.wn < 0.419) or (prop.wn > 420.0) then
      prop.kwx := 4;

    for j := 0 to 1 do
      if (prop.hg[j] < 0.5) or (prop.hg[j] > 3000.0) then
        prop.kwx := 4;

    static.dmin := abs(prop.he[0] - prop.he[1]) / 200e-3;
    q := adiff2(0.0, prop, propa, static);
    static.xae := power(prop.wn * (prop.gme * prop.gme), -THIRD);
    d3 := mymax(propa.dlsa, 1.3787 * static.xae + propa.dla);
    d4 := d3 + 2.7574 * static.xae;
    a3 := adiff2(d3, prop, propa, static);
    a4 := adiff2(d4, prop, propa, static);
    propa.emd := (a4 - a3) / (d4 - d3);
    propa.aed := a3 - propa.emd * d3;
  end;

  if (prop.mdp >= 0) then  (* if initializing the area mode *)
  begin
    prop.mdp := 0;  (* area mode is initialized *)
    prop.dist := d;
  end;

  if (prop.dist > 0.0) then
  begin
    if (prop.dist > 1000e3) then
      (* prop.dist being in meters, if greater than 1000 km, kwx=1 *)
      prop.kwx := mymax(prop.kwx, 1);

    if (prop.dist < static.dmin) then
      prop.kwx := mymax(prop.kwx, 3);

    if (prop.dist < 1e3) or (prop.dist > 2000e3) then
      prop.kwx := 4;
  end;

  if (prop.dist < propa.dlsa) then
  begin
    if (iw <= 0.0) then  (* if interval width is zero or less, used for area mode *)
    begin
      if (not static.wlos) then
      begin
        q := alos2(0.0, prop, propa);
        d2 := propa.dlsa;
        a2 := propa.aed + d2 * propa.emd;
        d0 := 1.908 * prop.wn * prop.he[0] * prop.he[1];

        if (propa.aed > 0.0) then
        begin
          prop.aref :=
            propa.aed + propa.emd * prop.dist;
        end
        else
        begin
          if (propa.aed = 0.0) then
          begin
            d0 := mymin(d0, 0.5 * propa.dla);
            d1 := d0 + 0.25 * (propa.dla - d0);
          end
          else   (* aed less than zero *)
          begin

            d1 := mymax(-propa.aed / propa.emd, 0.25 * propa.dla);
          end;
          a1 := alos2(d1, prop, propa);
          wq := False;

          if (d0 < d1) then
          begin
            a0 := alos2(d0, prop, propa);
            a2 := mymin(a2, alos2(d2, prop, propa));
            q := ln(d2 / d0);
            propa.ak2 :=
              mymax(0.0, ((d2 - d0) * (a1 - a0) -
              (d1 - d0) * (a2 - a0)) / ((d2 - d0) *
              ln(d1 / d0) - (d1 - d0) * q));
            wq := (propa.aed >= 0.0) or (propa.ak2 > 0.0);

            if (wq) then
            begin
              propa.ak1 :=
                (a2 - a0 - propa.ak2 * q) / (d2 - d0);

              if (propa.ak1 < 0.0) then
              begin
                propa.ak1 := 0.0;
                propa.ak2 :=
                  FORTRAN_DIM(a2, a0) / q;

                if (propa.ak2 = 0.0) then
                  propa.
                    ak1 :=
                    propa.emd;
              end;
            end;
          end;

          if (not wq) then
          begin
            propa.ak1 :=
              FORTRAN_DIM(a2, a1) / (d2 - d1);
            propa.ak2 := 0.0;

            if (propa.ak1 = 0.0) then
              propa.ak1 := propa.emd;
          end;
          propa.ael :=
            a2 - propa.ak1 * d2 - propa.ak2 * ln(d2);
          static.wlos := True;
        end;
      end;
    end
    else
    begin  (* for ITWOM point-to-point mode *)

      if (not static.wlos) then
      begin
        q := alos2(0.0, prop, propa);  (* coefficient setup *)
        static.wlos := True;
      end;

      if (prop.los = 1) then  (* if line of sight *)
      begin
        prop.aref := alos2(pd1, prop, propa);
      end
      else
      begin
        if (trunc(prop.dist - prop.dl[0]) = 0) then  (* if at 1st horiz *)
        begin
          prop.aref :=
            5.8 + alos2(pd1, prop, propa);
        end
        else
        if (trunc(prop.dist - prop.dl[0]) > 0.0) then
          (* if past 1st horiz *)
        begin
          q := adiff2(0.0, prop, propa, static);
          prop.aref := adiff2(pd1, prop, propa, static);
        end
        else
        begin
          prop.aref := 1.0;
        end;

      end;
    end;
  end;

  (* los and diff. range coefficents done. Starting troposcatter *)
  if (prop.dist <= 0.0) or (prop.dist >= propa.dlsa) then
  begin
    if (iw = 0.0) then  (* area mode *)
    begin
      if (not static.wscat) then
      begin
        q := ascat(0.0, prop, propa, static);
        d5 := propa.dla + 200e3;
        d6 := d5 + 200e3;
        a6 := ascat(d6, prop, propa, static);
        a5 := ascat(d5, prop, propa, static);

        if (a5 < 1000.0) then
        begin
          propa.ems := (a6 - a5) / 200e3;
          propa.dx :=
            mymax(propa.dlsa, mymax(propa.dla + 0.3 *
            static.xae * ln(47.7 * prop.wn),
            (a5 - propa.aed - propa.ems * d5) / (propa.emd - propa.ems)));

          propa.aes :=
            (propa.emd - propa.ems) * propa.dx + propa.aed;
        end
        else
        begin
          propa.ems := propa.emd;
          propa.aes := propa.aed;
          propa.dx := 10000000;
        end;
        static.wscat := True;
      end;

      if (prop.dist > propa.dx) then
      begin
        prop.aref := propa.aes + propa.ems * prop.dist;
      end
      else
      begin
        prop.aref := propa.aed + propa.emd * prop.dist;
      end;
    end
    else
    begin  (* ITWOM mode  q used to preset coefficients with zero input *)

      if (not static.wscat) then
      begin
        d5 := 0.0;
        d6 := 0.0;
        q := ascat(0.0, prop, propa, static);
        a6 := ascat(pd1, prop, propa, static);
        q := adiff2(0.0, prop, propa, static);
        a5 := adiff2(pd1, prop, propa, static);

        if (a5 <= a6) then
        begin
          propa.dx := 10000000;
          prop.aref := a5;
        end
        else
        begin
          propa.dx := propa.dlsa;
          prop.aref := a6;
        end;
        static.wscat := True;
      end;
    end;
  end;
  prop.aref := mymax(prop.aref, 0.0);
end;

function curve(const c1, c2, x1, x2, x3, de: double): double;
var
  temp1, temp2: double;
begin
  (* result := (c1+c2/(1.0+power((de-x2)/x3,2.0)))*power(de/x1,2.0)/(1.0+power(de/x1,2.0)); *)
  temp1 := (de - x2) / x3;
  temp2 := de / x1;

  temp1 := temp1 * temp1;
  temp2 := temp2 * temp2;

  Result := (c1 + c2 / (1.0 + temp1)) * temp2 / (1.0 + temp2);
end;

function avar(const zzt, zzl, zzc: double; var prop: prop_type; var propv: propv_type;
  var static: static_type): double;
const
  bv1: array[0..6] of double = (-9.67, -0.62, 1.26, -9.21, -0.62, -0.39, 3.15);
  bv2: array[0..6] of double = (12.7, 9.19, 15.5, 9.05, 9.19, 2.86, 857.9);
  xv1: array[0..6] of double =
    (144.9e3, 228.9e3, 262.6e3, 84.1e3, 228.9e3, 141.7e3, 2222.e3);
  xv2: array[0..6] of double =
    (190.3e3, 205.2e3, 185.2e3, 101.1e3, 205.2e3, 315.9e3, 164.8e3);
  xv3: array[0..6] of double =
    (133.8e3, 143.6e3, 99.8e3, 98.6e3, 143.6e3, 167.4e3, 116.3e3);
  bsm1: array[0..6] of double = (2.13, 2.66, 6.11, 1.98, 2.68, 6.86, 8.51);
  bsm2: array[0..6] of double = (159.5, 7.67, 6.65, 13.11, 7.16, 10.38, 169.8);
  xsm1: array[0..6] of double =
    (762.2e3, 100.4e3, 138.2e3, 139.1e3, 93.7e3, 187.8e3, 609.8e3);
  xsm2: array[0..6] of double =
    (123.6e3, 172.5e3, 242.2e3, 132.7e3, 186.8e3, 169.6e3, 119.9e3);
  xsm3: array[0..6] of double =
    (94.5e3, 136.4e3, 178.6e3, 193.5e3, 133.5e3, 108.9e3, 106.6e3);
  bsp1: array[0..6] of double = (2.11, 6.87, 10.08, 3.68, 4.75, 8.58, 8.43);
  bsp2: array[0..6] of double = (102.3, 15.53, 9.60, 159.3, 8.12, 13.97, 8.19);
  xsp1: array[0..6] of double =
    (636.9e3, 138.7e3, 165.3e3, 464.4e3, 93.2e3, 216.0e3, 136.2e3);
  xsp2: array[0..6] of double =
    (134.8e3, 143.7e3, 225.7e3, 93.1e3, 135.9e3, 152.0e3, 188.5e3);
  xsp3: array[0..6] of double =
    (95.6e3, 98.6e3, 129.7e3, 94.2e3, 113.4e3, 122.7e3, 122.9e3);
  bsd1: array[0..6] of double = (1.224, 0.801, 1.380, 1.000, 1.224, 1.518, 1.518);
  bzd1: array[0..6] of double = (1.282, 2.161, 1.282, 20., 1.282, 1.282, 1.282);
  bfm1: array[0..6] of double = (1.0, 1.0, 1.0, 1.0, 0.92, 1.0, 1.0);
  bfm2: array[0..6] of double = (0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
  bfm3: array[0..6] of double = (0.0, 0.0, 0.0, 0.0, 1.77, 0.0, 0.0);
  bfp1: array[0..6] of double = (1.0, 0.93, 1.0, 0.93, 0.93, 1.0, 1.0);
  bfp2: array[0..6] of double = (0.0, 0.31, 0.0, 0.19, 0.31, 0.0, 0.0);
  bfp3: array[0..6] of double = (0.0, 2.00, 0.0, 1.79, 2.00, 0.0, 0.0);

var
  rt, rl, avarv, q, vs, zt, zl, zc: double;
  sgt, yr, temp1, temp2: double;
  temp_klim: integer;
begin
  rt := 7.8;
  rl := 24.0;
  temp_klim := propv.klim - 1;

  if (propv.lvar > 0) then
  begin
    if (propv.klim <= 0) or (propv.klim > 7) then
    begin
      propv.klim := 5;
      temp_klim := 4;
      prop.kwx := mymax(prop.kwx, 2);
    end;
    with static do
    begin
      cv1 := bv1[temp_klim];
      cv2 := bv2[temp_klim];
      yv1 := xv1[temp_klim];
      yv2 := xv2[temp_klim];
      yv3 := xv3[temp_klim];
      csm1 := bsm1[temp_klim];
      csm2 := bsm2[temp_klim];
      ysm1 := xsm1[temp_klim];
      ysm2 := xsm2[temp_klim];
      ysm3 := xsm3[temp_klim];
      csp1 := bsp1[temp_klim];
      csp2 := bsp2[temp_klim];
      ysp1 := xsp1[temp_klim];
      ysp2 := xsp2[temp_klim];
      ysp3 := xsp3[temp_klim];
      csd1 := bsd1[temp_klim];
      zd := bzd1[temp_klim];
      cfm1 := bfm1[temp_klim];
      cfm2 := bfm2[temp_klim];
      cfm3 := bfm3[temp_klim];
      cfp1 := bfp1[temp_klim];
      cfp2 := bfp2[temp_klim];
      cfp3 := bfp3[temp_klim];
    end;

    if propv.lvar = 4 then
    begin
      static.kdv := propv.mdvar;
      static.ws := static.kdv >= 20;

      if (static.ws) then
        static.kdv := static.kdv - 20;

      static.w1 := static.kdv >= 10;

      if (static.w1) then
        static.kdv := static.kdv - 10;

      if (static.kdv < 0) or (static.kdv > 3) then
      begin
        static.kdv := 0;
        prop.kwx := mymax(prop.kwx, 2);
      end;
    end;

    if propv.lvar in [4, 3] then
    begin
      q := ln(0.133 * prop.wn);

      (* gm:=cfm1+cfm2/(power(cfm3*q,2.0)+1.0); *)
      (* gp:=cfp1+cfp2/(power(cfp3*q,2.0)+1.0); *)

      static.gm := static.cfm1 + static.cfm2 / ((static.cfm3 * q * static.cfm3 * q) + 1.0);
      static.gp := static.cfp1 + static.cfp2 / ((static.cfp3 * q * static.cfp3 * q) + 1.0);
    end;

    if propv.lvar in [4, 3, 2] then
    begin
      static.dexa :=
        sqrt(18e6 * prop.he[0]) + sqrt(18e6 * prop.he[1]) + power(
        (575.7e12 / prop.wn), THIRD);
    end;

    if propv.lvar in [4, 3, 2, 1] then
    begin
      if (prop.dist < static.dexa) then
        static.de := 130e3 * prop.dist / static.dexa
      else
        static.de := 130e3 + prop.dist - static.dexa;
    end;

    with static do
    begin
      vmd := curve(cv1, cv2, yv1, yv2, yv3, de);
      sgtm := curve(csm1, csm2, ysm1, ysm2, ysm3, de) * gm;
      sgtp := curve(csp1, csp2, ysp1, ysp2, ysp3, de) * gp;
      sgtd := sgtp * csd1;
      tgtd := (sgtp - sgtd) * zd;
    end;

    if (static.w1) then
      static.sgl := 0.0
    else
    begin
      q := (1.0 - 0.8 * exp(-prop.dist / 50e3)) * prop.dh * prop.wn;
      static.sgl := 10.0 * q / (q + 13.0);
    end;

    if (static.ws) then
      static.vs0 := 0.0
    else
    begin
      (* vs0=power(5.0+3.0*exp(-de/100e3),2.0); *)
      temp1 := (5.0 + 3.0 * exp(-static.de / 100e3));
      static.vs0 := temp1 * temp1;

    end;

    propv.lvar := 0;
  end;

  zt := zzt;
  zl := zzl;
  zc := zzc;

  case static.kdv of
    0:
    begin
      zt := zc;
      zl := zc;
    end;
    1:
    begin
      zl := zc;
    end;

    2:
    begin
      zl := zt;
    end;
  end;

  if (abs(zt) > 3.1) or (abs(zl) > 3.1) or (abs(zc) > 3.1) then
    prop.kwx := mymax(prop.kwx, 1);

  if (zt < 0.0) then
    sgt := static.sgtm
  else if (zt <= static.zd) then
    sgt := static.sgtp
  else
    sgt := static.sgtd + static.tgtd / zt;

  (* vs=vs0+power(sgt*zt,2.0)/(rt+zc*zc)+power(sgl*zl,2.0)/(rl+zc*zc); *)

  temp1 := sgt * zt;
  temp2 := static.sgl * zl;

  vs := static.vs0 + (temp1 * temp1) / (rt + zc * zc) + (temp2 * temp2) /
    (rl + zc * zc);

  if (static.kdv = 0) then
  begin
    yr := 0.0;
    propv.sgc := sqrt(sgt * sgt + static.sgl * static.sgl + vs);
  end
  else if (static.kdv = 1) then
  begin
    yr := sgt * zt;
    propv.sgc := sqrt(static.sgl * static.sgl + vs);
  end
  else if (static.kdv = 2) then
  begin
    yr := sqrt(sgt * sgt + static.sgl * static.sgl) * zt;
    propv.sgc := sqrt(vs);
  end
  else
  begin
    yr := sgt * zt + static.sgl * zl;
    propv.sgc := sqrt(vs);
  end;

  avarv := prop.aref - static.vmd - yr - propv.sgc * zc;

  if (avarv < 0.0) then
    avarv := avarv * (29.0 - avarv) / (29.0 - 10.0 * avarv);

  Result := avarv;
end;

procedure hzns(const pfl: array of double; var prop: prop_type);
(* Used only with ITM 1.2.2 *)
var
  wq: boolean;
  np: integer;
  xi, za, zb, qc, q, sb, sa: double;
  i: integer;
begin
  np := trunc(pfl[0]);
  xi := pfl[1];
  za := pfl[2] + prop.hg[0];
  zb := pfl[np + 2] + prop.hg[1];
  qc := 0.5 * prop.gme;
  q := qc * prop.dist;
  prop.the[1] := (zb - za) / prop.dist;
  prop.the[0] := prop.the[1] - q;
  prop.the[1] := -prop.the[1] - q;
  prop.dl[0] := prop.dist;
  prop.dl[1] := prop.dist;

  if (np >= 2) then
  begin
    sa := 0.0;
    sb := prop.dist;
    wq := True;

    for i := 1 to np - 1 do
    begin
      sa := sa + xi;
      sb := sb - xi;
      q := pfl[i + 2] - (qc * sa + prop.the[0]) * sa - za;

      if (q > 0.0) then
      begin
        prop.the[0] := prop.the[0] + (q / sa);
        prop.dl[0] := sa;
        wq := False;
      end;

      if (not wq) then
      begin
        q := pfl[i + 2] - (qc * sb + prop.the[1]) * sb - zb;

        if (q > 0.0) then
        begin
          prop.the[1] := prop.the[1] + (q / sb);
          prop.dl[1] := sb;
        end;
      end;
    end;
  end;
end;

procedure hzns2(const pfl: array of double; var prop: prop_type;
  const {%H-}propa: propa_type);
var
  wq: boolean;
  np, rp, i, j: integer;
  xi, za, zb, qc, q, sb, sa, dr, dshh: double;

begin
  np := trunc(pfl[0]);
  xi := pfl[1];
  za := pfl[2] + prop.hg[0];
  zb := pfl[np + 2] + prop.hg[1];
  prop.tiw := xi;
  prop.ght := za;
  prop.ghr := zb;
  qc := 0.5 * prop.gme;
  q := qc * prop.dist;
  prop.the[1] := arctan((zb - za) / prop.dist);
  prop.the[0] := (prop.the[1]) - q;
  prop.the[1] := -prop.the[1] - q;
  prop.dl[0] := prop.dist;
  prop.dl[1] := prop.dist;
  prop.hht := 0.0;
  prop.hhr := 0.0;
  prop.los := 1;

  if (np >= 2) then
  begin
    sa := 0.0;
    sb := prop.dist;
    wq := True;

    for j := 1 to np - 1 do
    begin
      sa := sa + xi;
      q := pfl[j + 2] - (qc * sa + prop.the[0]) * sa - za;

      if (q > 0.0) then
      begin
        prop.los := 0;
        prop.the[0] := prop.the[0] + (q / sa);
        prop.dl[0] := sa;
        prop.the[0] := mymin(prop.the[0], 1.569);
        prop.hht := pfl[j + 2];
        wq := False;
      end;
    end;

    if (not wq) then
    begin
      for i := 1 to np - 1 do
      begin
        sb := sb - xi;
        q := pfl[np + 2 - i] - (qc * (prop.dist - sb) + prop.the[1]) *
          (prop.dist - sb) - zb;
        if (q > 0.0) then
        begin
          prop.the[1] := prop.the[1] + (q / (prop.dist - sb));
          prop.the[1] := mymin(prop.the[1], 1.57);
          prop.the[1] :=
            mymax(prop.the[1], -1.568);
          prop.hhr := pfl[np + 2 - i];
          prop.dl[1] := mymax(0.0, prop.dist - sb);
        end;
      end;
      prop.the[0] :=
        arctan((prop.hht - za) / prop.dl[0]) - 0.5 * prop.gme * prop.dl[0];
      prop.the[1] :=
        arctan((prop.hhr - zb) / prop.dl[1]) - 0.5 * prop.gme * prop.dl[1];
    end;
  end;

  if ((prop.dl[1]) < (prop.dist)) then
  begin
    dshh := prop.dist - prop.dl[0] - prop.dl[1];

    if (trunc(dshh) = 0) then  (* one obstacle *)
    begin
      dr := prop.dl[1] / (1 + zb / prop.hht);
    end
    else   (* two obstacles *)
    begin

      dr := prop.dl[1] / (1 + zb / prop.hhr);
    end;
  end
  else     (* line of sight  *)
  begin

    dr := (prop.dist) / (1 + zb / za);
  end;
  rp := 2 + trunc(floor(0.5 + dr / xi));
  prop.rpl := rp;
  prop.rph := pfl[rp];
end;

procedure z1sq1(const z: array of double; const x1, x2: double; out z0, zn: double);
(* Used only with ITM 1.2.2 *)
var
  xn, xa, xb, x, a, b: double;
  i, n, ja, jb: integer;
begin
  xn := z[0];
  xa := trunc(FORTRAN_DIM(x1 / z[1], 0.0));
  xb := xn - trunc(FORTRAN_DIM(xn, x2 / z[1]));

  if (xb <= xa) then
  begin
    xa := FORTRAN_DIM(xa, 1.0);
    xb := xn - FORTRAN_DIM(xn, xb + 1.0);
  end;

  ja := trunc(xa);
  jb := trunc(xb);
  n := jb - ja;
  xa := xb - xa;
  x := -0.5 * xa;
  xb := xb + x;

  a := 0.5 * (z[ja + 2] + z[jb + 2]);
  b := 0.5 * (z[ja + 2] - z[jb + 2]) * x;

  for i := 2 to n do
  begin
    Inc(ja);
    x := x + 1.0;
    a := a + z[ja + 2];
    b := b + z[ja + 2] * x;
  end;

  a := a / xa;
  b := b * 12.0 / ((xa * xa + 2.0) * xa);
  z0 := a - b * xb;
  zn := a + b * (xn - xb);
end;

procedure z1sq2(const z: array of double; const x1, x2: double; out z0, zn: double);
(* corrected for use with ITWOM *)
var
  xn, xa, xb, x, a, b, bn: double;
  i, n, ja, jb: integer;
begin
  xn := z[0];
  xa := trunc(FORTRAN_DIM(x1 / z[1], 0.0));
  xb := xn - trunc(FORTRAN_DIM(xn, x2 / z[1]));

  if (xb <= xa) then
  begin
    xa := FORTRAN_DIM(xa, 1.0);
    xb := xn - FORTRAN_DIM(xn, xb + 1.0);
  end;

  ja := trunc(xa);
  jb := trunc(xb);
  xa := (2 * int((xb - xa) / 2)) - 1;
  x := -0.5 * (xa + 1);
  xb := xb + x;
  ja := int64(jb) - 1 - trunc(xa);
  n := jb - ja;
  a := (z[ja + 2] + z[jb + 2]);
  b := (z[ja + 2] - z[jb + 2]) * x;
  bn := 2 * (x * x);

  for i := 2 to n do
  begin
    Inc(ja);
    x := x + 1.0;
    bn := bn + (x * x);
    a := a + z[ja + 2];
    b := b + z[ja + 2] * x;
  end;

  a := a / (xa + 2);
  b := b / bn;
  z0 := a - (b * xb);
  zn := a + (b * (xn - xb));
end;

function qtile(const nn: integer; var _a: array of double; const offset: integer;
  const ir: integer): double;
var
  q, r: double;
  m, n, i, j, j1, i0, k: integer;
  done, goto10: boolean;
begin
  q := 0.0;  (* q initialization -- KD2BD *)
  j1 := 0;
  i0 := 0;  (* more initializations -- KD2BD *)
  done := False;
  goto10 := True;

  m := 0;
  n := nn;
  k := mymin(mymax(0, ir), n);

  while (not done) do
  begin
    if (goto10) then
    begin
      q := _a[offset + k];
      i0 := m;
      j1 := n;
    end;

    i := i0;

    while (i <= n) and (_a[offset + i] >= q) do
      Inc(i);

    if (i > n) then
      i := n;

    j := j1;

    while (j >= m) and (_a[offset + j] <= q) do
      Dec(j);

    if (j < m) then
      j := m;

    if (i < j) then
    begin
      r := _a[offset + i];
      _a[offset + i] := _a[offset + j];
      _a[offset + j] := r;
      i0 := i + 1;
      j1 := j - 1;
      goto10 := False;
    end
    else if (i < k) then
    begin
      _a[offset + k] := _a[offset + i];
      _a[offset + i] := q;
      m := i + 1;
      goto10 := True;
    end
    else if (j > k) then
    begin
      _a[offset + k] := _a[offset + j];
      _a[offset + j] := q;
      n := j - 1;
      goto10 := True;
    end
    else
      done := True;
  end;

  Result := q;
end;

function qerf(const z: double): double;

const
  b1 = 0.319381530;
  b2 = -0.356563782;
  b3 = 1.781477937;
  b4 = -1.821255987;
  b5 = 1.330274429;
  rp = 4.317008;
  rrt2pi = 0.398942280;
var
  t, x, qerfv: double;
begin
  x := z;
  t := abs(x);

  if (t >= 10.0) then
    qerfv := 0.0
  else
  begin
    t := rp / (t + rp);
    qerfv :=
      exp(-0.5 * x * x) * rrt2pi * ((((b5 * t + b4) * t + b3) * t + b2) *
      t + b1) * t;
  end;

  if (x < 0.0) then
    qerfv := 1.0 - qerfv;

  Result := qerfv;
end;

function d1thx(const pfl: array of double; const x1, x2: double): double;
var
  np, ka, kb, n, k, j: integer;
  d1thxv, sn, xa, xb: double;
  s: array of double;
begin
  np := trunc(pfl[0]);
  xa := x1 / pfl[1];
  xb := x2 / pfl[1];
  d1thxv := 0.0;

  if (xb - xa < 2.0) then  // exit out
  begin
    Result := d1thxv;
    exit;
  end;

  ka := trunc(0.1 * (xb - xa + 8.0));
  ka := mymin(mymax(4, ka), 25);
  n := 10 * ka - 5;
  kb := n - ka + 1;
  sn := n - 1;
  setlength(s, n + 2);
  s[0] := sn;
  s[1] := 1.0;
  xb := (xb - xa) / sn;
  k := trunc(xa + 1.0);
  xa := xa - k;

  for j := 0 to n - 1 do
  begin
    while (xa > 0.0) and (k < np) do
    begin
      xa := xa - 1.0;
      Inc(k);
    end;

    s[j + 2] := pfl[k + 2] + (pfl[k + 2] - pfl[k + 1]) * xa;
    xa := xa + xb;
  end;

  z1sq1(s, 0.0, sn, xa, xb);
  xb := (xb - xa) / sn;

  for j := 0 to n - 1 do
  begin
    s[j + 2] := s[j + 2] - xa;
    xa := xa + xb;
  end;

  d1thxv := qtile(n - 1, s, 2, ka - 1) - qtile(n - 1, s, 2, kb - 1);
  d1thxv := d1thxv / (1.0 - 0.8 * exp(-(x2 - x1) / 50.0e3));
  SetLength(s, 0);

  Result := d1thxv;
end;

function d1thx2(const pfl: array of double; const x1, x2: double;
  const {%H-}propa: propa_type): double;
var
  np, ka, kb, n, k, kmx, j: integer;
  d1thx2v, sn, xa, xb, xc: double;
  s: array of double;
begin
  np := trunc(pfl[0]);
  xa := x1 / pfl[1];
  xb := x2 / pfl[1];
  d1thx2v := 0.0;

  if (xb - xa < 2.0) then  // exit out
  begin
    Result := d1thx2v;
    exit;
  end;

  ka := trunc(0.1 * (xb - xa + 8.0));
  kmx := mymax(25, trunc(83350 / (pfl[1])));
  ka := mymin(mymax(4, ka), kmx);
  n := 10 * ka - 5;
  kb := n - ka + 1;
  sn := n - 1;
  setlength(s, n + 2);
  s[0] := sn;
  s[1] := 1.0;
  xb := (xb - xa) / sn;
  k := trunc(xa + 1.0);
  xc := xa - k;

  for j := 0 to n - 1 do
  begin
    while (xc > 0.0) and (k < np) do
    begin
      xc := xc - 1.0;
      Inc(k);
    end;

    s[j + 2] := pfl[k + 2] + (pfl[k + 2] - pfl[k + 1]) * xc;
    xc := xc + xb;
  end;

  z1sq2(s, 0.0, sn, xa, xb);
  xb := (xb - xa) / sn;

  for j := 0 to n - 1 do
  begin
    s[j + 2] := s[j + 2] - xa;
    xa := xa + xb;
  end;

  d1thx2v := qtile(n - 1, s, 2, ka - 1) - qtile(n - 1, s, 2, kb - 1);
  d1thx2v := d1thx2v / (1.0 - 0.8 * exp(-(x2 - x1) / 50.0e3));
  setlength(s, 0);
  Result := d1thx2v;
end;

procedure qlrpfl(const pfl: array of double; const klimx, mdvarx: integer;
  var prop: prop_type; var propa: propa_type; var propv: propv_type;
  var static: static_type);
var
  np, j: integer;
  xl: array[0..1] of double;
  q, za, zb, temp: double;
begin

  prop.dist := pfl[0] * pfl[1];
  np := trunc(pfl[0]);
  hzns(pfl, prop);

  for j := 0 to 1 do
    xl[j] := mymin(15.0 * prop.hg[j], 0.1 * prop.dl[j]);

  xl[1] := prop.dist - xl[1];
  prop.dh := d1thx(pfl, xl[0], xl[1]);

  if (prop.dl[0] + prop.dl[1] > 1.5 * prop.dist) then
  begin
    z1sq1(pfl, xl[0], xl[1], za, zb);
    prop.he[0] := prop.hg[0] + FORTRAN_DIM(pfl[2], za);
    prop.he[1] := prop.hg[1] + FORTRAN_DIM(pfl[np + 2], zb);

    for j := 0 to 1 do
      prop.dl[j] :=
        sqrt(2.0 * prop.he[j] / prop.gme) *
        exp(-0.07 * sqrt(prop.dh /
        mymax(prop.he[j], 5.0)));

    q := prop.dl[0] + prop.dl[1];

    if (q <= prop.dist) then
      (* if there is a rounded horizon, or two obstructions, in the path *)
    begin
      (* q=power(prop.dist/q,2.0); *)
      temp := prop.dist / q;
      q := temp * temp;

      for j := 0 to 1 do
      begin
        prop.he[j] := prop.he[j] * q;
        (* tx effective height set to be path dist/distance between obstacles *)
        prop.dl[j] :=
          sqrt(2.0 * prop.he[j] / prop.gme) * exp(-0.07 *
          sqrt(prop.dh / mymax(prop.he[j], 5.0)));
      end;
    end;

    for j := 0 to 1 do
      (* original empirical adjustment?  uses delta-h to adjust grazing angles *)
    begin
      q := sqrt(2.0 * prop.he[j] / prop.gme);
      prop.the[j] :=
        (0.65 * prop.dh * (q / prop.dl[j] - 1.0) - 2.0 * prop.he[j]) / q;
    end;
  end
  else
  begin
    z1sq1(pfl, xl[0], 0.9 * prop.dl[0], za, q);
    z1sq1(pfl, prop.dist - 0.9 * prop.dl[1], xl[1], q, zb);
    prop.he[0] := prop.hg[0] + FORTRAN_DIM(pfl[2], za);
    prop.he[1] := prop.hg[1] + FORTRAN_DIM(pfl[np + 2], zb);
  end;

  prop.mdp := -1;
  propv.lvar := mymax(propv.lvar, 3);

  if (mdvarx >= 0) then
  begin
    propv.mdvar := mdvarx;
    propv.lvar := mymax(propv.lvar, 4);
  end;

  if (klimx > 0) then
  begin
    propv.klim := klimx;
    propv.lvar := 5;
  end;

  lrprop(0.0, prop, propa, static);
end;

procedure qlrpfl2(const pfl: array of double; const klimx, mdvarx: integer;
  var prop: prop_type; var propa: propa_type; var propv: propv_type;
  var static: static_type);
var
  np, j: integer;
  xl: array[0..1] of double;
  dlb, q, za, zb, temp, rad, rae1, rae2: double;
begin
  prop.dist := pfl[0] * pfl[1];
  np := trunc(pfl[0]);
  hzns2(pfl, prop, propa);
  dlb := prop.dl[0] + prop.dl[1];
  prop.rch[0] := prop.hg[0] + pfl[2];
  prop.rch[1] := prop.hg[1] + pfl[np + 2];

  for j := 0 to 1 do
    xl[j] := mymin(15.0 * prop.hg[j], 0.1 * prop.dl[j]);

  xl[1] := prop.dist - xl[1];
  prop.dh := d1thx2(pfl, xl[0], xl[1], propa);

  if ((np < 1) or (pfl[1] > 150.0)) then
  begin
    (* for TRANSHORIZON; diffraction over a mutual horizon, or for one or more obstructions *)
    if (dlb < 1.5 * prop.dist) then
    begin
      z1sq2(pfl, xl[0], 0.9 * prop.dl[0], za, q);
      z1sq2(pfl, prop.dist - 0.9 * prop.dl[1], xl[1], q, zb);
      prop.he[0] := prop.hg[0] + FORTRAN_DIM(pfl[2], za);
      prop.he[1] := prop.hg[1] + FORTRAN_DIM(pfl[np + 2], zb);
    end
    (* for a Line-of-Sight path *)
    else
    begin
      z1sq2(pfl, xl[0], xl[1], za, zb);
      prop.he[0] := prop.hg[0] + FORTRAN_DIM(pfl[2], za);
      prop.he[1] := prop.hg[1] + FORTRAN_DIM(pfl[np + 2], zb);

      for j := 0 to 1 do
        prop.dl[j] :=
          sqrt(2.0 * prop.he[j] / prop.gme) * exp(-0.07 *
          sqrt(prop.dh / mymax(prop.he[j], 5.0)));

      (* for one or more obstructions only NOTE buried as in ITM FORTRAN and DLL, not functional  *)
      if ((prop.dl[0] + prop.dl[1]) <= prop.dist) then
      begin
        (* q=power(prop.dist/(dl[0]+dl[1])),2.0); *)
        temp := prop.dist / (prop.dl[0] + prop.dl[1]);
        q := temp * temp;
      end;

      for j := 0 to 1 do
      begin
        prop.he[j] := prop.he[j] * q;
        prop.dl[j] :=
          sqrt(2.0 * prop.he[j] / prop.gme) * exp(-0.07 *
          sqrt(prop.dh / mymax(prop.he[j], 5.0)));
      end;

      (* this sets (or resets) prop.the, and is not in The Guide FORTRAN QLRPFL *)
      for j := 0 to 1 do
      begin
        q := sqrt(2.0 * prop.he[j] / prop.gme);
        prop.the[j] :=
          (0.65 * prop.dh * (q / prop.dl[j] - 1.0) - 2.0 * prop.he[j]) / q;
      end;
    end;
  end
  else       (* for ITWOM ,computes he for tx, rcvr, and the receiver approach angles for use in saalos *)
  begin
    prop.he[0] := prop.hg[0] + (pfl[2]);
    prop.he[1] := prop.hg[1] + (pfl[np + 2]);

    rad := (prop.dist - 500.0);

    if (prop.dist > 550.0) then
    begin
      z1sq2(pfl, rad, prop.dist, rae1, rae2);
    end
    else
    begin
      rae1 := 0.0;
      rae2 := 0.0;
    end;

    prop.thera := arctan(abs(rae2 - rae1) / prop.dist);

    if (rae2 < rae1) then
    begin
      prop.thera := -prop.thera;
    end;

    prop.thenr :=
      arctan(mymax(0.0, (pfl[np + 2] - pfl[np + 1])) / pfl[1]);
  end;

  prop.mdp := -1;
  propv.lvar := mymax(propv.lvar, 3);

  if (mdvarx >= 0) then
  begin
    propv.mdvar := mdvarx;
    propv.lvar := mymax(propv.lvar, 4);
  end;

  if (klimx > 0) then
  begin
    propv.klim := klimx;
    propv.lvar := 5;
  end;

  lrprop2(0.0, prop, propa, static);
end;

function deg2rad(const d: double): double;
begin
  Result := d * 3.1415926535897 / 180.0;
end;

(***************************************************************************************
 * Point-To-Point Mode Calculations
 ***************************************************************************************)

procedure point_to_point_ITM(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref, frq_mhz: double;
  const radio_climate, pol: integer; const conf, rel: double;
  out dbloss: double; out strmode: string; out errnum: integer);

(******************************************************************************

Note that point_to_point has become point_to_point_ITM for use as the old ITM

  pol:
    0-Horizontal, 1-Vertical

  radio_climate:
    1-Equatorial, 2-Continental Subtropical,
    3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
    6-Maritime Temperate, Over Land, 7-Maritime Temperate,
    Over Sea

  conf, rel: .01 to .99

  elev[]: [num points - 1], [delta dist(meters)],
          [height(meters) point 1], ..., [height(meters) point n]

  errnum: 0- No Error.
    1- Warning: Some parameters are nearly out of range.
                Results should be used with caution.
    2- Note: Default parameters have been substituted for
             impossible ones.
    3- Warning: A combination of parameters is out of range.
          Results are probably invalid.
    Other-  Warning: Some parameters are out of range.
      Results are probably invalid.

*****************************************************************************)
var
  prop: prop_type;
  propv: propv_type;
  propa: propa_type;
  global: static_type;
  zsys, zc, zr: double;
  eno, enso, q: double;
  ja, jb, i, np: integer;
  (* dkm, xkm: double; *)
  fs: double;
begin
  fillchar({%H-}prop, sizeof(prop), 0);
  fillchar({%H-}propv, sizeof(propv), 0);
  fillchar({%H-}propa, sizeof(propa), 0);
  fillchar({%H-}global, sizeof(global), 0);

  zsys := 0;

  prop.hg[0] := tht_m;
  prop.hg[1] := rht_m;
  propv.klim := radio_climate;
  prop.kwx := 0;
  propv.lvar := 5;
  prop.mdp := -1;
  zc := qerfi(conf);
  zr := qerfi(rel);
  np := trunc(elev[0]);
  (* dkm:=(elev[1]*elev[0])/1000.0; *)
  (* xkm:=elev[1]/1000.0; *)
  eno := eno_ns_surfref;
  enso := 0.0;
  q := enso;

  if (q <= 0.0) then
  begin
    ja := trunc(3.0 + 0.1 * elev[0]);  (* added (long) to correct *)
    jb := np - ja + 6;

    for i := ja - 1 to jb - 1 do
      zsys := zsys + elev[i];

    zsys := zsys / (jb - ja + 1);
    q := eno;
  end;

  propv.mdvar := 12;
  qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
  qlrpfl(elev, propv.klim, propv.mdvar, prop, propa, propv, global);
  fs := 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
  q := prop.dist - propa.dla;

  if (trunc(q) < 0.0) then
  begin
    strmode := 'Line-Of-Sight Mode';
  end
  else
  begin
    if (trunc(q) = 0.0) then
      strmode := 'Single Horizon'
    else if (trunc(q) > 0.0) then
      strmode := 'Double Horizon';

    if (prop.dist <= propa.dlsa) or (prop.dist <= propa.dx) then
      strmode := strmode + ', Diffraction Dominant'
    else if (prop.dist > propa.dx) then
      strmode := strmode + ', Troposcatter Dominant';
  end;

  dbloss := avar(zr, 0.0, zc, prop, propv, global) + fs;
  errnum := prop.kwx;
end;

procedure point_to_point(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref, frq_mhz: double;
  const radio_climate, pol: integer; const conf, rel: double;
  out dbloss: double; out strmode: string; out errnum: integer);
(******************************************************************************

  Note that point_to_point_two has become point_to_point
  for drop-in interface to splat.cpp.
  The new variable inputs,
  double enc_ncc_clcref,
  double clutter_height,
  double clutter_density,
  double delta_h_diff, and
  int mode_var)
  have been given fixed values below.

  pol:
    0-Horizontal, 1-Vertical, 2-Circular

  radio_climate:
    1-Equatorial, 2-Continental Subtropical,
    3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
    6-Maritime Temperate, Over Land, 7-Maritime Temperate,
    Over Sea

  conf, rel: .01 to .99

  elev[]: [num points - 1], [delta dist(meters)],
          [height(meters) point 1], ..., [height(meters) point n]

  clutter_height    25.2 meters for compatibility with ITU-R P.1546-2.

  clutter_density   1.0 for compatibility with ITU-R P.1546-2.

  delta_h_diff    optional delta h for beyond line of sight. 90 m. average.
        setting to 0.0 will default to use of original internal
        use of delta-h for beyond line-of-sight range.

  mode_var    set to 12; or to 1 for FCC ILLR;  see documentation

  enc_ncc_clcref     clutter refractivity; 1000 N-units to match ITU-R P.1546-2

  eno=eno_ns_surfref  atmospheric refractivity at sea level; 301 N-units nominal
        (ranges from 250 for dry, hot day to 450 on hot, humid day]
        (stabilizes near 301 in cold, clear weather)

  errnum: 0- No Error.
    1- Warning: Some parameters are nearly out of range.
                Results should be used with caution.
    2- Note: Default parameters have been substituted for
             impossible ones.
    3- Warning: A combination of parameters is out of range.
          Results are probably invalid.
    Other-  Warning: Some parameters are out of range.
      Results are probably invalid.

*****************************************************************************)
var
  prop: prop_type;
  propv: propv_type;
  propa: propa_type;
  global: static_type;
  zsys, zc, zr: double;
  eno, enso, q: double;
  ja, jb, i, np: integer;
  (* dkm, xkm: double; *)
  tpd, fs: double;
  mode_var: integer;
begin
  fillchar({%H-}prop, sizeof(prop), 0);
  fillchar({%H-}propv, sizeof(propv), 0);
  fillchar({%H-}propa, sizeof(propa), 0);
  fillchar({%H-}global, sizeof(global), 0);

  zsys := 0;

  prop.hg[0] := tht_m;
  prop.hg[1] := rht_m;
  propv.klim := radio_climate;
  prop.kwx := 0;
  propv.lvar := 5;
  prop.mdp := -1;
  prop.ptx := pol;
  prop.thera := 0.0;
  prop.thenr := 0.0;
  zc := qerfi(conf);
  zr := qerfi(rel);
  np := trunc(elev[0]);
  (* dkm:=(elev[1]*elev[0])/1000.0; *)
  (* xkm:=elev[1]/1000.0; *)
  eno := eno_ns_surfref;
  enso := 0.0;
  q := enso;

  (* PRESET VALUES for Basic Version w/o additional inputs active *)

  prop.encc := 1000.00;  (*  double enc_ncc_clcref preset  *)
  prop.cch := 22.5;  (* double clutter_height preset to ILLR calibration.;
           use 25.3 for ITU-P1546-2 calibration *)
  prop.cd := 1.00;    (* double clutter_density preset *)
  mode_var := 1;  (* int mode_var set to 1 for FCC compatibility;
           normally, SPLAT presets this to 12 *)
  prop.dhd := 0.0;    (* delta_h_diff preset *)

  if (q <= 0.0) then
  begin
    ja := trunc(3.0 + 0.1 * elev[0]);
    jb := np - ja + 6;

    for i := ja - 1 to jb - 1 do
      zsys := zsys + elev[i];

    zsys := zsys / (jb - ja + 1);
    q := eno;
  end;

  propv.mdvar := mode_var;
  qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
  qlrpfl2(elev, propv.klim, propv.mdvar, prop, propa, propv, global);
  tpd :=
    sqrt((prop.he[0] - prop.he[1]) * (prop.he[0] - prop.he[1]) +
    (prop.dist) * (prop.dist));
  fs := 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(tpd / 1000.0);
  q := prop.dist - propa.dla;

  if (trunc(q) < 0.0) then
    strmode := 'L-o-S'
  else
  begin
    if (trunc(q) = 0.0) then
      strmode := '1_Hrzn'
    else if (trunc(q) > 0.0) then
      strmode := '2_Hrzn';

    if (prop.dist <= propa.dlsa) or (prop.dist <= propa.dx) then
    begin
      if trunc(prop.dl[1]) = 0.0 then
        strmode := strmode + '_Peak'
      else
        strmode := strmode + '_Diff';

    end
    else if (prop.dist > propa.dx) then
      strmode := strmode + '_Tropo';
  end;

  dbloss := avar(zr, 0.0, zc, prop, propv, global) + fs;
  errnum := prop.kwx;
end;

procedure point_to_pointMDH_two(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
  enc_ncc_clcref, clutter_height, clutter_density, delta_h_diff, frq_mhz: double;
  const radio_climate, pol: integer; const timepct, locpct, confpct: double;
  out dbloss: double; out propmode: integer; out deltaH: double; out errnum: integer);
(*************************************************************************************************
   pol: 0-Horizontal, 1-Vertical
   radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
                  4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
                  7-Maritime Temperate, Over Sea
   timepct, locpct, confpct: .01 to .99
   elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
   propmode:  Value   Mode
               -1     mode is undefined
                0     Line of Sight
                5     Single Horizon, Diffraction
                6     Single Horizon, Troposcatter
                9     Double Horizon, Diffraction
               10     Double Horizon, Troposcatter
   errnum: 0- No Error.
           1- Warning: Some parameters are nearly out of range.
                       Results should be used with caution.
           2- Note: Default parameters have been substituted for impossible ones.
           3- Warning: A combination of parameters is out of range.
                       Results are probably invalid.
           Other-  Warning: Some parameters are out of range.
                            Results are probably invalid.
*************************************************************************************************)
var
  prop: prop_type;
  propv: propv_type;
  propa: propa_type;
  global: static_type;
  zsys: double;
  eno, enso, q: double;
  ja, jb, i, np: integer;
  (* dkm, xkm: double; *)
  fs: double;

  ztime, zloc, zconf: double;
begin
  fillchar({%H-}prop, sizeof(prop), 0);
  fillchar({%H-}propv, sizeof(propv), 0);
  fillchar({%H-}propa, sizeof(propa), 0);
  fillchar({%H-}global, sizeof(global), 0);

  zsys := 0;

  propmode := -1;    // mode is undefined
  prop.hg[0] := tht_m;
  prop.hg[1] := rht_m;
  propv.klim := radio_climate;
  prop.encc := enc_ncc_clcref;
  prop.cch := clutter_height;
  prop.cd := clutter_density;
  prop.dhd := delta_h_diff;
  prop.kwx := 0;
  propv.lvar := 5;
  prop.mdp := -1;
  prop.ptx := pol;
  prop.thera := 0.0;
  prop.thenr := 0.0;
  ztime := qerfi(timepct);
  zloc := qerfi(locpct);
  zconf := qerfi(confpct);
  np := trunc(elev[0]);
  (* dkm := (elev[1] * elev[0]) / 1000.0; *)
  (* xkm := elev[1] / 1000.0; *)
  eno := eno_ns_surfref;
  enso := 0.0;
  q := enso;

  (* PRESET VALUES for Basic Version w/o additional inputs active *)

  prop.encc := 1000.00;  (*  double enc_ncc_clcref  *)
  prop.cch := 22.5;  (* double clutter_height *)
  prop.cd := 1.00;    (* double clutter_density *)

  if (q <= 0.0) then
  begin
    ja := trunc(3.0 + 0.1 * elev[0]);  (* to match addition of (long) *)
    jb := np - ja + 6;
    for i := ja - 1 to jb - 1 do
      zsys := zsys + elev[i];
    zsys := zsys / (jb - ja + 1);
    q := eno;
  end;
  propv.mdvar := 12;
  qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
  qlrpfl2(elev, propv.klim, propv.mdvar, prop, propa, propv, global);
  fs := 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);

  deltaH := prop.dh;
  q := prop.dist - propa.dla;
  if (trunc(q) < 0.0) then
    propmode := 0  // L-of-S
  else
  begin
    if (trunc(q) = 0.0) then
      propmode := 4  // 1-Hrzn
    else if (trunc(q) > 0.0) then
      propmode := 8;  // 2-Hrzn
    if (prop.dist <= propa.dlsa) or (prop.dist <= propa.dx) then
      propmode := propmode + 1  // Diff
    else if (prop.dist > propa.dx) then
      propmode := propmode + 2;  // Tropo
  end;
  dbloss := avar(ztime, zloc, zconf, prop, propv, global) + fs;
  //avar(time,location,confidence)
  errnum := prop.kwx;
end;

procedure point_to_pointDH(const elev: array of double;
  const tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref,
  enc_ncc_clcref, clutter_height, clutter_density, delta_h_diff, frq_mhz: double;
  const radio_climate, pol: integer; const conf, rel: double;
  out dbloss, deltaH: double; out errnum: integer);
(*************************************************************************************************
   pol: 0-Horizontal, 1-Vertical
   radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
                  4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
                  7-Maritime Temperate, Over Sea
   conf, rel: .01 to .99
   elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
   errnum: 0- No Error.
           1- Warning: Some parameters are nearly out of range.
                       Results should be used with caution.
           2- Note: Default parameters have been substituted for impossible ones.
           3- Warning: A combination of parameters is out of range.
                       Results are probably invalid.
           Other-  Warning: Some parameters are out of range.
                            Results are probably invalid.
*************************************************************************************************)
var
  strmode: string;
  prop: prop_type;
  propv: propv_type;
  propa: propa_type;
  global: static_type;
  zc, zr: double;
  eno, enso, q: double;
  ja, jb, i, np: longint;
  (* dkm, xkm: double; *)
  fs: double;
  zsys: double;
begin
  fillchar({%H-}prop, sizeof(prop), 0);
  fillchar({%H-}propv, sizeof(propv), 0);
  fillchar({%H-}propa, sizeof(propa), 0);
  fillchar({%H-}global, sizeof(global), 0);
  zsys := 0;

  prop.hg[0] := tht_m;
  prop.hg[1] := rht_m;
  propv.klim := radio_climate;
  prop.encc := enc_ncc_clcref;
  prop.cch := clutter_height;
  prop.cd := clutter_density;
  prop.dhd := delta_h_diff;
  prop.kwx := 0;
  propv.lvar := 5;
  prop.mdp := -1;
  prop.ptx := pol;
  prop.thera := 0.0;
  prop.thenr := 0.0;
  zc := qerfi(conf);
  zr := qerfi(rel);
  np := trunc(elev[0]);
  (* dkm := (elev[1] * elev[0]) / 1000.0; *)
  (* xkm := elev[1] / 1000.0; *)
  eno := eno_ns_surfref;
  enso := 0.0;
  q := enso;

  (* PRESET VALUES for Basic Version w/o additional inputs active *)

  prop.encc := 1000.00;  (*  double enc_ncc_clcref  *)
  prop.cch := 22.5;  (* double clutter_height *)
  prop.cd := 1.00;    (* double clutter_density *)

  if (q <= 0.0) then
  begin
    ja := trunc(3.0 + 0.1 * elev[0]);  (* to match KD2BD addition of (long)  *)
    jb := np - ja + 6;
    for i := ja - 1 to jb - 1 do
      zsys := zsys + elev[i];
    zsys := zsys / (jb - ja + 1);
    q := eno;
  end;
  propv.mdvar := 12;
  qlrps(frq_mhz, zsys, q, pol, eps_dielect, sgm_conductivity, prop);
  qlrpfl2(elev, propv.klim, propv.mdvar, prop, propa, propv, global);
  fs := 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
  deltaH := prop.dh;
  q := prop.dist - propa.dla;
  if (trunc(q) < 0.0) then
    strmode := 'Line-Of-Sight Mode'
  else
  begin
    if (trunc(q) = 0.0) then
      strmode := 'Single Horizon'
    else if (trunc(q) > 0.0) then
      strmode := 'Double Horizon';
    if (prop.dist <= propa.dlsa) or (prop.dist <= propa.dx) then
      strmode := strmode + ', Diffraction Dominant'
    else if (prop.dist > propa.dx) then
      strmode := strmode + ', Troposcatter Dominant';
  end;
  dbloss := avar(zr, 0.0, zc, prop, propv, global) + fs;  //avar(time,location,confidence)
  errnum := prop.kwx;
end;

(********************************************************
 * Area Mode Calculations                               *
 ********************************************************)
procedure area(const ModVar: longint; const deltaH, tht_m, rht_m, dist_km: double;
  const TSiteCriteria, RSiteCriteria: integer;
  const eps_dielect, sgm_conductivity, eno_ns_surfref, enc_ncc_clcref,
  clutter_height, clutter_density, delta_h_diff, frq_mhz: double;
  const radio_climate, pol: integer; const pctTime, pctLoc, pctConf: double;
  out dbloss: double; out errnum: integer);

// pol: 0-Horizontal, 1-Vertical
// TSiteCriteria, RSiteCriteria:
//                 0 - random, 1 - careful, 2 - very careful

// radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
//                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
//                7-Maritime Temperate, Over Sea
// ModVar: 0 - Single: pctConf is "Time/Situation/Location", pctTime, pctLoc not used
//         1 - Individual: pctTime is "Situation/Location", pctConf is "Confidence", pctLoc not used
//         2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is "Confidence", pctLoc not used
//         3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is "Confidence"
// pctTime, pctLoc, pctConf: .01 to .99
// errnum: 0- No Error.
//         1- Warning: Some parameters are nearly out of range.
//                     Results should be used with caution.
//         2- Note: Default parameters have been substituted for impossible ones.
//         3- Warning: A combination of parameters is out of range.
//                     Results are probably invalid.
//         Other-  Warning: Some parameters are out of range.
//                          Results are probably invalid.
// NOTE: strmode is not used at this time.
var
  prop: prop_type;
  propv: propv_type;
  propa: propa_type;
  global: static_type;
  zt, zl, zc, xlb: double;
  fs: double;
  ivar: integer;
  eps, eno, sgm: double;
  ipol: integer;
  kst: array[0..1] of integer;

begin
  fillchar({%H-}prop, sizeof(prop), 0);
  fillchar({%H-}propv, sizeof(propv), 0);
  fillchar({%H-}propa, sizeof(propa), 0);
  fillchar({%H-}global, sizeof(global), 0);

  kst[0] := TSiteCriteria;
  kst[1] := RSiteCriteria;
  zt := qerfi(pctTime / 100.0);
  zl := qerfi(pctLoc / 100.0);
  zc := qerfi(pctConf / 100.0);
  eps := eps_dielect;
  sgm := sgm_conductivity;
  eno := eno_ns_surfref;
  prop.dh := deltaH;
  prop.hg[0] := tht_m;
  prop.hg[1] := rht_m;
  propv.klim := radio_climate;
  prop.encc := enc_ncc_clcref;
  prop.cch := clutter_height;
  prop.cd := clutter_density;
  prop.dhd := delta_h_diff;
  prop.ens := eno;
  prop.kwx := 0;
  ivar := ModVar;
  ipol := pol;
  qlrps(frq_mhz, 0.0, eno, ipol, eps, sgm, prop);
  qlra(kst, propv.klim, ivar, prop, propv);

  if (propv.lvar < 1) then
    propv.lvar := 1;

  lrprop2(dist_km * 1000.0, prop, propa, global);
  fs := 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
  xlb := fs + avar(zt, zl, zc, prop, propv, global);
  dbloss := xlb;
  if (prop.kwx = 0) then
    errnum := 0
  else
    errnum := prop.kwx;
end;

function ITMAreadBLoss(const ModVar: integer;
  const deltaH, tht_m, rht_m, dist_km: double;
  const TSiteCriteria, RSiteCriteria: integer;
  const eps_dielect, sgm_conductivity, eno_ns_surfref,
  enc_ncc_clcref, clutter_height, clutter_density, delta_h_diff,
  frq_mhz: double;
  const radio_climate, pol: integer;
  const pctTime, pctLoc, pctConf: double): double;
var
  errnum: integer;
  dbLoss: double;
begin
  area(ModVar, deltaH, tht_m, rht_m, dist_km,
    TSiteCriteria, RSiteCriteria,
    eps_dielect, sgm_conductivity, eno_ns_surfref, enc_ncc_clcref,
    clutter_height, clutter_density, delta_h_diff, frq_mhz,
    radio_climate, pol, pctTime, pctLoc, pctConf,
    dbloss, errnum);

  Result := dbloss;
end;

// other models

function _10log10f(const x: double): double; inline;
begin
  Result := 4.342944 * ln(x);
end;

function _20log10f(x: double): double; inline;
begin
  Result := 8.685889 * ln(x);
end;

function COST231pathLoss(f, TxH, RxH, d: double; const mode: TEnvironment): double;
(*
COST231 extension to HATA model
Frequency 1500 to 2000MHz
TxH = Base station height 30 to 200m
RxH = Mobile station height 1 to 10m
Distance 1-20km
modes 1 = URBAN, 2 = SUBURBAN, 3 = OPEN
http://morse.colorado.edu/~tlen5510/text/classwebch3.html
*)
var
  C: integer;
  lRxH, C_H: double;
  c0, cf: double;
  logf: double;
begin
  C := 3;    // 3dB for Urban
  lRxH := log10(11.75 * RxH);
  C_H := 3.2 * (lRxH * lRxH) - 4.97;  // Large city (conservative)
  c0 := 69.55;
  cf := 26.16;
  if (f > 1500) then
  begin
    c0 := 46.3;
    cf := 33.9;
  end;
  if (mode = TEnvironment.evSuburban) then
  begin
    C := 0;    // Medium city (average)
    lRxH := log10(1.54 * RxH);
    C_H := 8.29 * (lRxH * lRxH) - 1.1;
  end
  else
  if (mode = TEnvironment.evUrban) then
  begin
    C := -3;    // Small city (Optimistic)
    C_H := (1.1 * log10(f) - 0.7) * RxH - (1.56 * log10(f)) + 0.8;
  end;
  logf := log10(f);
  Result :=
    c0 + (cf * logf) - (13.82 * log10(TxH)) - C_H + (44.9 - 6.55 * log10(TxH)) *
    log10(d) + C;
end;

function ECC33pathLoss(f, TxH, RxH, d: double; const mode: TEnvironment): double;
var
  Gr, Gb: double;
  Afs, Abm: double;
begin

  // Sanity check as this model operates within limited Txh/Rxh bounds
  if (TxH - RxH < 0) then
    RxH := RxH / (d * 2);

  // MHz to GHz
  f := f / 1000;

  Gr := 0.759 * RxH - 1.862;  // Big city with tall buildings (1)
  // PL := Afs + Abm - Gb - Gr
  Afs := 92.4 + 20 * log10(d) + 20 * log10(f);
  Abm :=
    20.41 + 9.83 * log10(d) + 7.894 * log10(f) + 9.56 * (log10(f) * log10(f));
  Gb := log10(TxH / 200) * (13.958 + 5.8 * (log10(d) * log10(d)));
  if (mode in [TEnvironment.evSuburban, TEnvironment.evRural]) then
    // Medium city (Europe)
    Gr := (42.57 + 13.7 * log10(f)) * (log10(RxH) - 0.585);

  Result := Afs + Abm - Gb - Gr;
end;


function EgliPathLoss(f, h1, h2, d: double): double;
var
  Lp50: double;
  C1, C2: double;
begin
  if (h1 > 10.0) and (h2 > 10.0) then
  begin
    Lp50 := 85.9;
    C1 := 2.0;
    C2 := 2.0;
  end
  else if (h1 > 10.0) then
  begin
    Lp50 := 76.3;
    C1 := 2.0;
    C2 := 1.0;
  end
  else if (h2 > 10.0) then
  begin
    Lp50 := 76.3;
    C1 := 1.0;
    C2 := 2.0;
  end
  else  // both antenna heights below 10 metres
  begin
    Lp50 := 66.7;
    C1 := 1.0;
    C2 := 1.0;
  end;  // end if

  Result := Lp50 + (4.0 * _10log10f(d) + 2.0 * _10log10f(f) - C1 *
    _10log10f(h1) - C2 * _10log10f(h2));
end;

function EricssonpathLoss(f, TxH, RxH, d: double; const Mode: TEnvironment): double;
var
  a0, a1, a2, a3: double;
  g1, g2: double;
begin
  (*
     AKA Ericsson 9999 model
   *)
  // Urban
  a0 := 36.2;
  a1 := 30.2;
  a2 := -12;
  a3 := 0.1;
  if (mode = TEnvironment.evSuburban) then
  begin
    // Suburban / Med loss
    a0 := 43.2;
    a1 := 68.93;
  end
  else
  if (mode = TEnvironment.evRural) then
  begin
    // Rural
    a0 := 45.95;
    a1 := 100.6;
  end;
  g1 := 3.2 * (log10(11.75 * RxH) * log10(11.75 * RxH));
  g2 := 44.49 * log10(f) - 4.78 * (log10(f) * log10(f));

  Result := a0 + a1 * log10(d) + a2 * log10(TxH) + a3 * log10(TxH) * log10(d) - g1 + g2;
end;

function FSPLpathLoss(f, d: double): double;
begin
  Result := (32.44 + _20log10f(f) + _20log10f(d));
end;

function HATApathLoss(f, h_B, h_M, d: double; mode: TEnvironment): double;
var
  L_u, lh_M, C_H: double;
  logf: double;
  logf_28: double;
begin
(*
HATA URBAN model for cellular planning
Frequency (MHz) 150 to 1500MHz
Base station height 30-200m
Mobile station height 1-10m
Distance 1-20km
mode 1 = URBAN
mode 2 = SUBURBAN
mode 3 = OPEN
*)
  logf := log10(f);

  if (f < 200) then
  begin
    lh_M := log10(1.54 * h_M);
    C_H := 8.29 * (lh_M * lh_M) - 1.1;
  end
  else
  begin
    lh_M := log10(11.75 * h_M);
    C_H := 3.2 * (lh_M * lh_M) - 4.97;
  end;

  L_u := 69.55 + 26.16 * logf - 13.82 * log10(h_B) - C_H +
    (44.9 - 6.55 * log10(h_B)) * log10(d);

  case mode of
    TEnvironment.evUrban:
    begin // Urban
      Result := L_u;
    end;
    TEnvironment.evSuburban:
    begin  // SubUrban
      logf_28 := log10(f / 28);
      Result := L_u - 2 * logf_28 * logf_28 - 5.4;
      exit;
    end;
    TEnvironment.evRural:
    begin // Open
      Result := L_u - 4.78 * logf * logf + 18.33 * logf - 40.94;
    end;
    else
      Result := 0;
  end;
end;

function PlaneEarthLoss(d, TxH, RxH: double): double;
begin
(*
Plane Earth Loss model
Frequency: N/A
Distance (km): Any
*)
  // Plane earth loss is independent of frequency.
  Result := 40 * log10(d) + 20 * log10(TxH) + 20 * log10(RxH);
end;

function SoilPathLoss(f, d, terdic: double): double;
var
  soil: double;
begin
  soil := (120 / terdic);
  Result := (6.4 + _20log10f(d) + _20log10f(f) + (8.69 * soil));
end;

function SUIpathLoss(f, TxH, RxH, d: double; mode: TEnvironment): double;
var
  a, b, c, s, XhCF: double;
  d0, _A, y: double;
  Xf, Xh: double;
begin
        (*
           f = Frequency (MHz) 1900 to 11000
           TxH =  Transmitter height (m)
           RxH = Receiver height (m)
           d = distance (km)
           mode A1 = URBAN / OBSTRUCTED
           mode B2 = SUBURBAN / PARTIALLY OBSTRUCTED
           mode C3 = RURAL / OPEN
           Paper 1 has a Rx height correction of / 2000
           Paper 2 has the same correction as / 2 and gives better results
           "Ranked number 2 University in the wurld"
           http://www.cl.cam.ac.uk/research/dtg/lce-pub/public/vsa23/VTC05_Empirical.pdf
           https://mentor.ieee.org/802.19/file/08/19-08-0010-00-0000-sui-path-loss-model.doc
         *)
  d := d * 1e3;               // km to m

  // Urban (A1) is default
  a := 4.6;
  b := 0.0075;
  c := 12.6;
  s := 8.2; // Optional fading value. 8.2 to 10.6dB
  XhCF := -10.8;

  if (mode = TEnvironment.evSuburban) then  // Suburban
  begin
    a := 4.0;
    b := 0.0065;
    c := 17.1;
    XhCF := -10.8;
  end;
  if (mode = TEnvironment.evRural) then  // Rural
  begin
    a := 3.6;
    b := 0.005;
    c := 20;
    XhCF := -20;
  end;
  d0 := 100.0;
  _A := _20log10f((4 * PI * d0) / (300.0 / f));
  y := a - (b * TxH) + (c / TxH);

  // Assume 2.4GHz
  Xf := 0;
  Xh := 0;

  //Correction factors for > 2GHz
  if (f > 2000) then
  begin
    Xf := 6.0 * log10(f / 2.0);
    Xh := XhCF * log10(RxH / 2.0);
  end;
  Result := _A + (10 * y) * (log10(d / d0)) + Xf + Xh + s;
end;

end.














































































