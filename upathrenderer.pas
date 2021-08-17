unit upathrenderer;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, upathloss, bgrabitmap, BGRABitmapTypes;

type
  TRegion = record
    color: array[0..127, 0..2] of integer;
    level: array[0..127] of integer;
    levels: integer;
  end;

  TColorRange = record
    dfVal: double;
    nR, nG, nB: byte;
  end;

  TRenderType = (rtSignal, rtPathLoss, rtRXPower, rtLineOfSight, rtElevation);

  TPathRenderer = class
  private
    FSignalRegion: TRegion;
    FLossRegion: TRegion;
    FDBmRegion: TRegion;
    FColorSet: array of TColorRange;

    procedure LoadSignalColors;
    function FindMatch(const APathLoss: TPathloss; const x, y: integer;
      const ARenderType: TRenderType; const AContourThreshold: integer;
      const ARegion: TRegion; out DrawContour: boolean): byte;

    function hillshade(const APathLoss: TPathLoss; const lat, lon: double): double;
    function Colorize(dfVal, dfHillshade: double): TBGRAPixel;
    procedure LoadDefaultColorSet;
  public
    constructor Create;
    function Render(const APathLoss: TPathloss; const ABitmap: TBGRABitmap;
      const ARenderType: TRenderType; const APlotTerrain: boolean = True;
      const AContourThreshold: integer = 0): boolean;
  end;

implementation

uses Math;

const
  GAMMA = 2.5;

constructor TPathRenderer.Create;
begin
  LoadSignalColors;
  LoadDefaultColorSet;
end;

procedure TPathRenderer.LoadSignalColors;
begin

  FSignalRegion.level[0] := 128;
  FSignalRegion.color[0][0] := 255;
  FSignalRegion.color[0][1] := 0;
  FSignalRegion.color[0][2] := 0;

  FSignalRegion.level[1] := 118;
  FSignalRegion.color[1][0] := 255;
  FSignalRegion.color[1][1] := 165;
  FSignalRegion.color[1][2] := 0;

  FSignalRegion.level[2] := 108;
  FSignalRegion.color[2][0] := 255;
  FSignalRegion.color[2][1] := 206;
  FSignalRegion.color[2][2] := 0;

  FSignalRegion.level[3] := 98;
  FSignalRegion.color[3][0] := 255;
  FSignalRegion.color[3][1] := 255;
  FSignalRegion.color[3][2] := 0;

  FSignalRegion.level[4] := 88;
  FSignalRegion.color[4][0] := 184;
  FSignalRegion.color[4][1] := 255;
  FSignalRegion.color[4][2] := 0;

  FSignalRegion.level[5] := 78;
  FSignalRegion.color[5][0] := 0;
  FSignalRegion.color[5][1] := 255;
  FSignalRegion.color[5][2] := 0;

  FSignalRegion.level[6] := 68;
  FSignalRegion.color[6][0] := 0;
  FSignalRegion.color[6][1] := 208;
  FSignalRegion.color[6][2] := 0;

  FSignalRegion.level[7] := 58;
  FSignalRegion.color[7][0] := 0;
  FSignalRegion.color[7][1] := 196;
  FSignalRegion.color[7][2] := 196;

  FSignalRegion.level[8] := 48;
  FSignalRegion.color[8][0] := 0;
  FSignalRegion.color[8][1] := 148;
  FSignalRegion.color[8][2] := 255;

  FSignalRegion.level[9] := 38;
  FSignalRegion.color[9][0] := 80;
  FSignalRegion.color[9][1] := 80;
  FSignalRegion.color[9][2] := 255;

  FSignalRegion.level[10] := 28;
  FSignalRegion.color[10][0] := 0;
  FSignalRegion.color[10][1] := 38;
  FSignalRegion.color[10][2] := 255;

  FSignalRegion.level[11] := 18;
  FSignalRegion.color[11][0] := 142;
  FSignalRegion.color[11][1] := 63;
  FSignalRegion.color[11][2] := 255;

  FSignalRegion.level[12] := 8;
  FSignalRegion.color[12][0] := 140;
  FSignalRegion.color[12][1] := 0;
  FSignalRegion.color[12][2] := 128;

  FSignalRegion.levels := 13;

  // PathLoss
  FLossRegion.level[0] := 80;
  FLossRegion.color[0][0] := 255;
  FLossRegion.color[0][1] := 0;
  FLossRegion.color[0][2] := 0;

  FLossRegion.level[1] := 90;
  FLossRegion.color[1][0] := 255;
  FLossRegion.color[1][1] := 128;
  FLossRegion.color[1][2] := 0;

  FLossRegion.level[2] := 100;
  FLossRegion.color[2][0] := 255;
  FLossRegion.color[2][1] := 165;
  FLossRegion.color[2][2] := 0;

  FLossRegion.level[3] := 110;
  FLossRegion.color[3][0] := 255;
  FLossRegion.color[3][1] := 206;
  FLossRegion.color[3][2] := 0;

  FLossRegion.level[4] := 120;
  FLossRegion.color[4][0] := 255;
  FLossRegion.color[4][1] := 255;
  FLossRegion.color[4][2] := 0;

  FLossRegion.level[5] := 130;
  FLossRegion.color[5][0] := 184;
  FLossRegion.color[5][1] := 255;
  FLossRegion.color[5][2] := 0;

  FLossRegion.level[6] := 140;
  FLossRegion.color[6][0] := 0;
  FLossRegion.color[6][1] := 255;
  FLossRegion.color[6][2] := 0;

  FLossRegion.level[7] := 150;
  FLossRegion.color[7][0] := 0;
  FLossRegion.color[7][1] := 208;
  FLossRegion.color[7][2] := 0;

  FLossRegion.level[8] := 160;
  FLossRegion.color[8][0] := 0;
  FLossRegion.color[8][1] := 196;
  FLossRegion.color[8][2] := 196;

  FLossRegion.level[9] := 170;
  FLossRegion.color[9][0] := 0;
  FLossRegion.color[9][1] := 148;
  FLossRegion.color[9][2] := 255;

  FLossRegion.level[10] := 180;
  FLossRegion.color[10][0] := 80;
  FLossRegion.color[10][1] := 80;
  FLossRegion.color[10][2] := 255;

  FLossRegion.level[11] := 190;
  FLossRegion.color[11][0] := 0;
  FLossRegion.color[11][1] := 38;
  FLossRegion.color[11][2] := 255;

  FLossRegion.level[12] := 200;
  FLossRegion.color[12][0] := 142;
  FLossRegion.color[12][1] := 63;
  FLossRegion.color[12][2] := 255;

  FLossRegion.level[13] := 210;
  FLossRegion.color[13][0] := 196;
  FLossRegion.color[13][1] := 54;
  FLossRegion.color[13][2] := 255;

  FLossRegion.level[14] := 220;
  FLossRegion.color[14][0] := 255;
  FLossRegion.color[14][1] := 0;
  FLossRegion.color[14][2] := 255;

  FLossRegion.level[15] := 230;
  FLossRegion.color[15][0] := 255;
  FLossRegion.color[15][1] := 194;
  FLossRegion.color[15][2] := 204;

  FLossRegion.levels := 16;

  FDBmRegion.level[0] := 0;
  FDBmRegion.color[0][0] := 255;
  FDBmRegion.color[0][1] := 0;
  FDBmRegion.color[0][2] := 0;

  FDBmRegion.level[1] := -10;
  FDBmRegion.color[1][0] := 255;
  FDBmRegion.color[1][1] := 128;
  FDBmRegion.color[1][2] := 0;

  FDBmRegion.level[2] := -20;
  FDBmRegion.color[2][0] := 255;
  FDBmRegion.color[2][1] := 165;
  FDBmRegion.color[2][2] := 0;

  FDBmRegion.level[3] := -30;
  FDBmRegion.color[3][0] := 255;
  FDBmRegion.color[3][1] := 206;
  FDBmRegion.color[3][2] := 0;

  FDBmRegion.level[4] := -40;
  FDBmRegion.color[4][0] := 255;
  FDBmRegion.color[4][1] := 255;
  FDBmRegion.color[4][2] := 0;

  FDBmRegion.level[5] := -50;
  FDBmRegion.color[5][0] := 184;
  FDBmRegion.color[5][1] := 255;
  FDBmRegion.color[5][2] := 0;

  FDBmRegion.level[6] := -60;
  FDBmRegion.color[6][0] := 0;
  FDBmRegion.color[6][1] := 255;
  FDBmRegion.color[6][2] := 0;

  FDBmRegion.level[7] := -70;
  FDBmRegion.color[7][0] := 0;
  FDBmRegion.color[7][1] := 208;
  FDBmRegion.color[7][2] := 0;

  FDBmRegion.level[8] := -80;
  FDBmRegion.color[8][0] := 0;
  FDBmRegion.color[8][1] := 196;
  FDBmRegion.color[8][2] := 196;

  FDBmRegion.level[9] := -90;
  FDBmRegion.color[9][0] := 0;
  FDBmRegion.color[9][1] := 148;
  FDBmRegion.color[9][2] := 255;

  FDBmRegion.level[10] := -100;
  FDBmRegion.color[10][0] := 80;
  FDBmRegion.color[10][1] := 80;
  FDBmRegion.color[10][2] := 255;

  FDBmRegion.level[11] := -110;
  FDBmRegion.color[11][0] := 0;
  FDBmRegion.color[11][1] := 38;
  FDBmRegion.color[11][2] := 255;

  FDBmRegion.level[12] := -120;
  FDBmRegion.color[12][0] := 142;
  FDBmRegion.color[12][1] := 63;
  FDBmRegion.color[12][2] := 255;

  FDBmRegion.level[13] := -130;
  FDBmRegion.color[13][0] := 196;
  FDBmRegion.color[13][1] := 54;
  FDBmRegion.color[13][2] := 255;

  FDBmRegion.level[14] := -140;
  FDBmRegion.color[14][0] := 255;
  FDBmRegion.color[14][1] := 0;
  FDBmRegion.color[14][2] := 255;

  FDBmRegion.level[15] := -150;
  FDBmRegion.color[15][0] := 255;
  FDBmRegion.color[15][1] := 194;
  FDBmRegion.color[15][2] := 204;

  FDBmRegion.levels := 16;
end;

procedure TPathRenderer.LoadDefaultColorSet;

  procedure Add(const AElevation: double; const AR, AG, AB: byte);
  begin
    SetLength(FColorSet, Length(FColorSet) + 1);
    with FColorSet[high(FColorSet)] do
    begin
      dfVal := AElevation;
      nR := AR;
      nG := AG;
      nB := AB;
    end;
  end;

begin
  SetLength(FColorSet, 0);
  Add(0, 102, 153, 153);
  Add(1, 46, 154, 88);
  Add(600, 251, 255, 128);
  Add(1200, 224, 108, 31);
  Add(2500, 200, 55, 55);
  Add(4000, 215, 244, 244);
end;



function TPathRenderer.FindMatch(const APathLoss: TPathloss;
  const x, y: integer; const ARenderType: TRenderType;
  const AContourThreshold: integer; const ARegion: TRegion;
  out DrawContour: boolean): byte;
var
  signal, loss, dBm: integer;
  Data: byte;
  z: integer;
begin
  Data := APathLoss.GetSignalValue(y * APathLoss.Height + x);

  case ARenderType of
    rtElevation:
    begin
      Result := 0;
      DrawContour := False;
      exit;
    end;
    rtSignal:
    begin
      signal := Data - 100;
      DrawContour := (AContourThreshold <> 0) and (signal < AContourThreshold);
      if (signal >= ARegion.level[0]) then
      begin
        Result := 0;
        exit;
      end
      else
      begin
        for z := 1 to ARegion.levels - 1 do
        begin
          if (signal < ARegion.level[z - 1]) and (signal >= ARegion.level[z]) then
          begin
            Result := z;
            exit;
          end;
        end;
      end;
    end;
    rtPathLoss:
    begin
      loss := Data;
      DrawContour := (Loss = 0) or ((AContourThreshold <> 0) and
        (loss > abs(AContourThreshold)));
      if (loss <= ARegion.level[0]) then
      begin
        Result := 0;
        exit;
      end
      else
      begin
        for z := 1 to ARegion.levels - 1 do
          if (z < ARegion.levels) then
            if (loss >= ARegion.level[z - 1]) and (loss < ARegion.level[z]) then
            begin
              Result := z;
              exit;
            end;
      end;
    end;
    rtRXPower:
    begin
      dBm := Data - 200;
      DrawContour := (AContourThreshold <> 0) and (dBm < AContourThreshold);
      if (dBm >= ARegion.level[0]) then
      begin
        Result := 0;
        exit;
      end
      else
      begin
        for z := 1 to ARegion.levels - 1 do
          if (z < ARegion.levels) then
            if (dBm < ARegion.level[z - 1]) and (dBm >= ARegion.level[z]) then
            begin
              Result := z;
              exit;
            end;
      end;
    end;
  end;
  DrawContour := False;
  Result := 255;
end;

function TPathRenderer.Colorize(dfVal, dfHillshade: double): TBGRAPixel;
var
  i, upper, mid, lower: integer;
  dfRatio: double;

  function Fix(d: double): byte;
  var
    i: integer;
  begin
    i := round(d);
    if i < 0 then
      i := 0;
    if i > 255 then
      i := 255;
    Result := byte(i);
  end;

begin
  lower := 0;
  dfVal := dfVal / 3.28084;
  Result.FromRGB(0, 0, 0);

  i := 0;
  upper := Length(FColorSet) - 1;
  while True do
  begin
    mid := (lower + upper) div 2;
    if (upper - lower <= 1) then
    begin
      if (dfVal <= FColorSet[lower].dfVal) then
        i := lower
      else if (dfVal <= FColorSet[upper].dfVal) then
        i := upper
      else
        i := upper + 1;
      break;
    end
    else if (FColorSet[mid].dfVal >= dfVal) then
      upper := mid
    else
      lower := mid;
  end;

  if i = 0 then
  begin
    Result.Red := FColorSet[0].nR;
    Result.Green := FColorSet[0].nG;
    Result.Blue := FColorSet[0].nB;
  end
  else
  begin
    if i >= Length(FColorSet) then
      i := Length(FColorSet) - 1;

    dfRatio :=
      (dfVal - FColorSet[i - 1].dfVal) / (FColorSet[i].dfVal -
      FColorSet[i - 1].dfVal);
    Result.Red := Fix(0.45 + FColorSet[i - 1].nR + dfRatio *
      (FColorSet[i].nR - FColorSet[i - 1].nR));
    Result.Green := Fix(0.45 + FColorSet[i - 1].nG + dfRatio *
      (FColorSet[i].nG - FColorSet[i - 1].nG));
    Result.Blue := Fix(0.45 + FColorSet[i - 1].nB + dfRatio *
      (FColorSet[i].nB - FColorSet[i - 1].nB));
  end;

  if (dfHillshade >= 0) and (dfHillshade <= 1) then
  begin
    Result.Red := round(dfHillshade * Result.Red);
    Result.Green := round(dfHillshade * Result.Green);
    Result.Blue := round(dfHillshade * Result.Blue);
  end;
end;

function TPathRenderer.hillshade(const APathLoss: TPathLoss;
  const lat, lon: double): double;

const
  Z_FACTOR = 1;
  KERNELSIZE = 1;
  ALTITUDE = 45 * PI / 180;
  AZIMUTH = 135 * PI / 180;
  ZENITH = PI / 2 - ALTITUDE;

var
  dpp, a, b, c, d, f, g, h, i: double;
  dzdx, dzdy, slope, aspect: double;
begin
  // Values in the eight neighboring cells
  with APathLoss do
  begin
    dpp := 1 / APathLoss.PixelPerDegree;
    GetElevation(self, lat - dpp, lon - dpp, a);
    GetElevation(self, lat, lon - dpp, b);
    GetElevation(self, lat + dpp, lon - dpp, c);
    GetElevation(self, lat - dpp, lon, d);
    GetElevation(self, lat + dpp, lon, f);
    GetElevation(self, lat - dpp, lon + dpp, g);
    GetElevation(self, lat, lon + dpp, h);
    GetElevation(self, lat + dpp, lon + dpp, i);
  end;
  dzdx := ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * KERNELSIZE);
  dzdy := ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * KERNELSIZE);
  slope := arctan(Z_FACTOR * sqrt(dzdx * dzdx + dzdy * dzdy));
  aspect := arctan2(dzdy, -dzdx);
  Result := ((cos(ZENITH) * cos(slope)) + (sin(ZENITH) * sin(slope) *
    cos(AZIMUTH - aspect)));
  if Result < 0 then
    Result := 0;
end;

function TPathRenderer.Render(const APathLoss: TPathloss;
  const ABitmap: TBGRABitmap; const ARenderType: TRenderType;
  const APlotTerrain: boolean = True; const AContourThreshold: integer = 0): boolean;

var
  red, green, blue: byte;
  mask: byte;

  offset: integer;
  x, y, x0, y0, match: integer;
  conversion, one_over_gamma: double;
  dpp: double;
  lat, lon: double;
  pelevation: double;
  drawcontour: boolean;
  plotterrain: boolean;
  Region: TRegion;

  procedure ADD_PIXEL(const x, y: integer; const r, g, b: byte);
  var
    c: TBGRAPixel;
  begin
    c.FromRGB(r, g, b);
    ABitmap.SetPixel(x, y, c);
  end;

begin
  Result := False;
  if not assigned(ABitmap) then
    exit;

  case ARenderType of
    rtSignal: Region := FSignalRegion;
    rtPathLoss: Region := FLossRegion;
    rtRXPower: Region := FDBmRegion;
    rtLineOfSight: Region := FSignalRegion;
  end;

  with APathLoss do
  begin
    if (MaxNorth=MinNorth) or (MaxWest=MinWest) then exit;
    drawcontour := False;
    dpp := 1 / PixelPerDegree;

    if (MaxElevation = -32768) or (MinElevation = 32768) then
      exit;

    ABitmap.SetSize(Height, Height);
    ABitmap.FillTransparent;

    one_over_gamma := 1.0 / GAMMA;
    conversion := Power((MaxElevation - MinElevation), one_over_gamma);
    if conversion > 0 then
      conversion := 255.0 / conversion;
    y := 0;
    repeat
      lat := MinNorth + (dpp * y);
      Inc(y);

      x := 0;
      repeat
        lon := MinWest + (dpp * x);
        Inc(x);

        if not SiteToOffset(lat, lon, x0, y0, offset) then
          continue;
        mask := GetMaskValue(offset);

        GetElevation(self, lat, lon, pelevation);
        red := 0;
        green := 0;
        blue := 0;

        if ARenderType = rtElevation then
        begin
          plotterrain := APlotTerrain;
          drawcontour := True;
        end
        else
        begin
          plotterrain := False;
          match := FindMatch(APathLoss, x0, y0, ARenderType,
            AContourThreshold, Region, drawcontour);

          if (match < Region.levels) then
          begin
            red := Region.color[match][0];
            green := Region.color[match][1];
            blue := Region.color[match][2];
          end;

          if (mask and 2 > 0) then
          begin
            (* Text Labels: Red or otherwise *)
            if (red >= 180) and (green <= 75) and (blue <= 75) then
              ADD_PIXEL(x0, y0, 255 xor red,
                255 xor green,
                255 xor blue)
            else
              ADD_PIXEL(x0, y0, 255, 0, 0);
          end
          else if (mask and 4 > 0) then
          begin
            (* County Boundaries: Black *)
            ADD_PIXEL(x0, y0, 0, 0, 0);
          end
          else
          begin
            if ARenderType = rtLineOfSight then
            begin
              case mask and 57 of
                1: ADD_PIXEL(x0, y0, 0, 255, 0); (* TX1: Green *)
                8: ADD_PIXEL(x0, y0, 0, 255, 255); (* TX2: Cyan *)
                9: ADD_PIXEL(x0, y0, 255, 255, 0);  (* TX1 + TX2: Yellow *)
                16: ADD_PIXEL(x0, y0, 147, 112, 219); (* TX3: Medium Violet *)
                17: ADD_PIXEL(x0, y0, 255, 192, 203); (* TX1 + TX3: Pink *)
                24: ADD_PIXEL(x0, y0, 255, 165, 0); (* TX2 + TX3: Orange *)
                25: ADD_PIXEL(x0, y0, 0, 100, 0); (* TX1 + TX2 + TX3: Dark Green *)
                32: ADD_PIXEL(x0, y0, 255, 130, 71); (* TX4: Sienna 1 *)
                33: ADD_PIXEL(x0, y0, 173, 255, 47); (* TX1 + TX4: Green Yellow *)
                40: ADD_PIXEL(x0, y0, 193, 255, 193); (* TX2 + TX4: Dark Sea Green 1 *)
                41: ADD_PIXEL(x0, y0, 255, 235, 205);
                (* TX1 + TX2 + TX4: Blanched Almond *)
                48: ADD_PIXEL(x0, y0, 0, 206, 209); (* TX3 + TX4: Dark Turquoise *)
                49: ADD_PIXEL(x0, y0, 0, 250, 154);
                (* TX1 + TX3 + TX4: Medium Spring Green *)
                56: ADD_PIXEL(x0, y0, 210, 180, 140); (* TX2 + TX3 + TX4: Tan *)
                57: ADD_PIXEL(x0, y0, 238, 201, 0); (* TX1 + TX2 + TX3 + TX4: Gold2 *)
                else
                  plotterrain := APlotTerrain;
              end;
            end
            else
            begin
              if (not drawcontour) and ((red <> 0) or (green <> 0) or (blue <> 0)) then
                ADD_PIXEL(x0, y0, red, green, blue)
              else
                plotterrain := APlotTerrain;
            end;
          end;
        end;
        if plotterrain then
          ABitmap.SetPixel(x0, y0,
            Colorize(pElevation, hillshade(APathLoss, lat, lon)));
      until lon > MaxWest;
    until lat > MaxNorth;
  end;
  Result := True;
end;

end.
