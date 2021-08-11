unit upathrenderer;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, upathloss, bgrabitmap;

type
  TRegion = record
    color: array[0..127, 0..2] of integer;
    level: array[0..127] of integer;
    levels: integer;
  end;

  TRenderType = (rtSignal, rtPathLoss, rtRXPower, rtLineOfSight);

  TPathRenderer = class
  private
    FSignalRegion: TRegion;
    FLossRegion: TRegion;
    FDBmRegion: TRegion;

    procedure LoadSignalColors;
    function FindMatch(const APathLoss: TPathloss; const x, y: integer;
      const ARenderType: TRenderType; const AContourThreshold: integer;
      const ARegion: TRegion;
      out DrawContour: boolean): byte;
  public
    constructor create;
    function Render(const APathLoss: TPathloss; const ABitmap: TBGRABitmap;
      const ARenderType: TRenderType; const APlotTerrain: boolean = True;
      const AContourThreshold: integer = 0): Boolean;
  end;

implementation

uses Math, BGRABitmapTypes;

const
  GAMMA = 2.5;

  constructor TPathRenderer.create;
  begin
    LoadSignalColors;
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

function TPathRenderer.FindMatch(const APathLoss: TPathloss;
  const x, y: integer; const ARenderType: TRenderType; const AContourThreshold: integer;
  const ARegion: TRegion;
  out DrawContour: boolean): byte;
var
  signal, loss, dBm: integer;
  data: byte;
  z: integer;
begin
    data := APathLoss.GetSignalValue(y * APathLoss.Height + x);

    case ARenderType of
      rtSignal:
      begin
        signal := data - 100;
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
        loss := data;
        DrawContour := (Loss=0) or ((AContourThreshold <> 0) and (loss > abs(AContourThreshold)));
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
        dBm := data - 200;
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

function TPathRenderer.Render(const APathLoss: TPathloss;
  const ABitmap: TBGRABitmap; const ARenderType: TRenderType;
  const APlotTerrain: boolean = True; const AContourThreshold: integer = 0): Boolean;

var
  terrain: byte;
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
  result := false;
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
    drawcontour := False;
    dpp := 1 / PixelPerDegree;

    if (MaxElevation = -32768) or (MinElevation = 32768) then exit;

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

          if plotterrain then
          begin
            if (pelevation = 0) then
              ADD_PIXEL(x0, y0, 0, 0, 170)
            else
            begin
              terrain :=
                byte(round(0.5  + power((pelevation - MinElevation), one_over_gamma) *
                conversion));
              ADD_PIXEL(x0, y0,
                terrain,
                terrain,
                terrain);
            end;
          end;
        end;
      until lon > MaxWest;
    until lat > MaxNorth;
  end;
  result := true;
end;

end.
