
from __future__ import annotations
import math
import csv
from dataclasses import dataclass
from typing import Dict, Tuple, List, Iterable
import os

A_EQ = 6378249.145
INV_F = 293.465
F = 1.0 / INV_F
E2 = F * (2 - F)
K0 = 0.9996
FALSE_EASTING = 500000.0
FALSE_NORTHING_S = 10000000.0

RAD = math.pi / 180.0
DEG = 180.0 / math.pi

TM_AUX_A = 6367386.644
TM_AUX_B = 16300.7
TM_AUX_C = 17.387
TM_AUX_E = 0.023

def dms_to_rad(d: int, m: int, s: float, sign: int = 1) -> float:
    vdeg = abs(d) + m/60.0 + s/3600.0
    return sign * vdeg * RAD

def rad_to_dms(theta: float):
    sign = 1 if theta >= 0 else -1
    v = abs(theta) * DEG
    d = int(math.floor(v + 1e-12))
    v = (v - d) * 60.0
    m = int(math.floor(v + 1e-12))
    s = (v - m) * 60.0
    if s >= 59.9999995:
        s = 0.0
        m += 1
    if m >= 60:
        m = 0
        d += 1
    return d, m, s, sign

def dms_str(theta: float, decimals: int = 2) -> str:
    d, m, s, sign = rad_to_dms(theta)
    sign_char = "-" if sign < 0 else ""
    return f"{sign_char}{d:02d}°{m:02d}'{s:0{3+decimals}.{decimals}f}\""

def to_polar(dx: float, dy: float):
    r = math.hypot(dx, dy)
    theta = math.atan2(dy, dx)
    return r, theta

def to_rect(r: float, theta: float):
    return r * math.cos(theta), r * math.sin(theta)

def bearing(ay: float, ax: float, by: float, bx: float) -> float:
    dx = bx - ax
    dy = by - ay
    _, th = to_polar(dx, dy)
    return th

def distance(ay: float, ax: float, by: float, bx: float) -> float:
    return math.hypot(bx - ax, by - ay)

def polygon_area(points):
    pts = list(points)
    if pts and pts[0] != pts[-1]:
        pts = pts + [pts[0]]
    acc = 0.0
    for (y1, x1), (y2, x2) in zip(pts, pts[1:]):
        acc += x1*y2 - x2*y1
    return 0.5 * acc

from dataclasses import dataclass

@dataclass
class Point:
    y: float
    x: float

class PointDB:
    def __init__(self, path: str):
        self.path = path
        self.points: Dict[str, Point] = {}
        if os.path.exists(path):
            self.load()
    def load(self):
        import csv
        with open(self.path, newline="", encoding="utf-8") as f:
            r = csv.reader(f)
            for row in r:
                if not row: continue
                pid, y, x = row
                self.points[pid] = Point(float(y), float(x))
    def save(self):
        import csv
        with open(self.path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            for pid, pt in self.points.items():
                w.writerow([pid, f"{pt.y:.4f}", f"{pt.x:.4f}"])
    def get(self, pid: str):
        return self.points.get(pid)
    def set(self, pid: str, y: float, x: float):
        self.points[pid] = Point(y, x)
    def delete(self, pid: str) -> bool:
        return self.points.pop(pid, None) is not None
    def prefix(self, pref: str):
        return {k: v for k, v in self.points.items() if k.startswith(pref)}

@dataclass
class HelmertParams:
    a: float
    b: float
    x0: float
    y0: float
    @property
    def scale(self) -> float:
        return math.hypot(self.a, self.b)
    @property
    def rotation(self) -> float:
        return math.atan2(self.b, self.a)

def fit_helmert(source, target) -> HelmertParams:
    if len(source) != len(target) or len(source) < 2:
        raise ValueError("Need at least 2 matching points with equal length.")
    n = len(source)
    ys = sum(y for y, _ in source) / n
    xs = sum(x for _, x in source) / n
    yt = sum(y for y, _ in target) / n
    xt = sum(x for _, x in target) / n
    src_c = [(y - ys, x - xs) for y, x in source]
    tgt_c = [(y - yt, x - xt) for y, x in target]
    Sxx = Syy = Sxy = Syx = 0.0
    denom = 0.0
    for (ysc, xsc), (ytc, xtc) in zip(src_c, tgt_c):
        Sxx += xsc * xtc
        Syy += ysc * ytc
        Sxy += xsc * ytc
        Syx += ysc * xtc
        denom += xsc*xsc + ysc*ysc
    a = (Sxx + Syy) / denom
    b = (Sxy - Syx) / denom
    x0 = xt - (a*xs - b*ys)
    y0 = yt - (b*xs + a*ys)
    return HelmertParams(a=a, b=b, x0=x0, y0=y0)

def apply_helmert(y: float, x: float, p: HelmertParams):
    x_p = p.a * x - p.b * y + p.x0
    y_p = p.b * x + p.a * y + p.y0
    return y_p, x_p

def slope_sea_level_correction(s: float, ha: float, hb: float, am: float = 6360000.0):
    c = -s * (ha + hb) / (2.0 * am) - ((ha - hb) ** 2) / (2.0 * s) + (s**3) / (25.6 * am**2)
    return c, s + c

def meridional_arc(phi: float) -> float:
    return (TM_AUX_A * phi
            - TM_AUX_B * math.sin(2*phi)
            + TM_AUX_C * math.sin(4*phi)
            - TM_AUX_E * math.sin(6*phi))

def footpoint(phi: float, M: float) -> float:
    for _ in range(20):
        f = TM_AUX_A*phi - TM_AUX_B*math.sin(2*phi) + TM_AUX_C*math.sin(4*phi) - TM_AUX_E*math.sin(6*phi) - M
        df = TM_AUX_A - 2*TM_AUX_B*math.cos(2*phi) + 4*TM_AUX_C*math.cos(4*phi) - 6*TM_AUX_E*math.cos(6*phi)
        step = f / df
        phi_new = phi - step
        if abs(step) < 1e-12:
            return phi_new
        phi = phi_new
    return phi

def utm_forward(lat: float, lon: float, lon0: float, south: bool=False):
    e2 = E2
    n2 = e2 / (1.0 - e2)
    def N(phi): return A_EQ / math.sqrt(1 - e2 * (math.sin(phi)**2))
    M = meridional_arc(lat)
    phi_f = footpoint(lat, M)
    T = math.tan(phi_f)**2
    C = n2 * (math.cos(phi_f)**2)
    A = (lon - lon0) * math.cos(phi_f)
    Nf = N(phi_f)
    E = FALSE_EASTING + K0 * Nf * (A + (1 - T + C) * (A**3) / 6.0 + (5 - 18*T + T*T + 72*C) * (A**5) / 120.0)
    Nn = K0 * (M + Nf * math.tan(phi_f) * ((A**2)/2.0 + (5 - T + 9*C + 4*C*C) * (A**4)/24.0
                                           + (61 - 58*T + T*T + 600*C) * (A**6)/720.0))
    if south:
        Nn += FALSE_NORTHING_S
    return (E, Nn)

def utm_inverse(E: float, Nn: float, lon0: float, south: bool=False):
    if south:
        Nn -= FALSE_NORTHING_S
    M = Nn / K0
    phi_f = footpoint(M / A_EQ, M / 1.0)
    e2 = E2
    n2 = e2 / (1.0 - e2)
    def N(phi): return A_EQ / math.sqrt(1 - e2 * (math.sin(phi)**2))
    Nf = N(phi_f)
    Rf = A_EQ * (1 - e2) / ((1 - e2 * (math.sin(phi_f)**2)) ** 1.5)
    D = (E - FALSE_EASTING) / (K0 * Nf)
    T = math.tan(phi_f)**2
    C = n2 * (math.cos(phi_f)**2)
    lat = (phi_f - (Nf * math.tan(phi_f) / Rf) *
           ( (D*D)/2.0 - (5 + 3*T + 10*C) * (D**4)/24.0
             + (61 + 90*T + 298*C + 45*T*T) * (D**6)/720.0 ))
    lon = lon0 + (D - (1 + 2*T + C) * (D**3)/6.0
                  + (5 + 28*T + 24*T*T + 6*C) * (D**5)/120.0) / math.cos(phi_f)
    return (lat, lon)

def compute_open_traverse(stations, distances, bearings):
    coords = [stations[0]]
    y, x = stations[0]
    for d, b in zip(distances, bearings):
        dx, dy = to_rect(d, b)
        x += dx
        y += dy
        coords.append((y, x))
    return coords

def bowditch_adjust(coords):
    n = len(coords)
    if n == 0:
        return coords
    if coords[0] != coords[-1]:
        coords = coords + [coords[0]]
        n += 1
    sum_len = 0.0
    for (y1, x1), (y2, x2) in zip(coords, coords[1:]):
        sum_len += math.hypot(x2-x1, y2-y1)
    mis_y = coords[-1][0] - coords[0][0]
    mis_x = coords[-1][1] - coords[0][1]
    adjusted = [coords[0]]
    acc_len = 0.0
    y, x = coords[0]
    for (y1, x1), (y2, x2) in zip(coords, coords[1:]):
        L = math.hypot(x2-x1, y2-y1)
        acc_len += L
        cy = -mis_y * (acc_len / sum_len)
        cx = -mis_x * (acc_len / sum_len)
        y = y2 + cy
        x = x2 + cx
        adjusted.append((y, x))
    return adjusted

#question the "cool shit"
def lo29_to_utm (E_lo: float,
                N_lo: float,
                *,
                fe_lo: float = 0.0,         # Lo false easting (if any)
                fn_lo: float = 0.0,         # Lo false northing (if any)
                k0_lo: float = 1.0,         # Lo central scale factor (often 1.0)
                south_hemisphere: bool = True,
                return_latlon: bool = False):
    """
    Convert a Lo29 Transverse Mercator coordinate (E_lo, N_lo) to UTM.

    Steps:
      1) Undo Lo false origins & scale, then invert TM @ 29°E to (lat, lon)
         using the same TM series as the Zambas port.
      2) Pick UTM zone from lon, then forward-project to UTM (k0=0.9996, FE=500000,
         FN=10,000,000 if south).

    Notes:
      - If your Lo system uses a false origin, set fe_lo/fn_lo accordingly.
      - If Lo uses a central scale k0 different from 1.0, set k0_lo.
      - Datum must match your Zambas math (this port assumes Intl 1924-style constants).
    """
    # --- 1) Prepare Lo coordinates for the Zambas TM-inverse call ---
    # utm_inverse expects: E_in includes its own 500,000 FE; N_in is scaled by K0 (0.9996).
    # To “trick” it into a generic TM inverse for Lo:
    #   - remove Lo FE/FN,
    #   - replace Lo FE with the UTM FE (add 500,000) so the internal subtraction matches,
    #   - rescale northing so M = N/k0_lo matches utm_inverse’s M = N_in/K0  ->  N_in = N * (K0/k0_lo)
    E_raw = E_lo - fe_lo
    N_raw = N_lo - fn_lo
    E_in  = E_raw + FALSE_EASTING                  # neutralize the built-in FE subtraction
    N_in  = N_raw * (K0 / k0_lo)                   # neutralize the built-in K0 division

    lon0_lo29 = 29.0 * RAD

    # IMPORTANT: for the Lo inverse, set south=False so utm_inverse does NOT subtract UTM's south FN.
    lat, lon = utm_inverse(E_in, N_in, lon0_lo29, south=False)

    # --- 2) Geographic -> UTM (proper) ---
    zone = int((lon * DEG + 180) // 6) + 1
    lon0_utm = (zone * 6 - 183) * RAD
    E_utm, N_utm = utm_forward(lat, lon, lon0_utm, south=south_hemisphere)

    if return_latlon:
        return {
            "lat_deg": lat * DEG,
            "lon_deg": lon * DEG,
            "utm_zone": zone,
            "E_utm": E_utm,
            "N_utm": N_utm
        }
    else:
        return zone, E_utm, N_utm
    

def lo27_to_utm27 (Y: float,X: float):
    """
    Converts Lo27 Transverse Mercator coordinates to UTM zone 27 coordinates.

    Parameters:
        E_lo (float): Easting in Lo27 coordinate system (meters).
        N_lo (float): Northing in Lo27 coordinate system (meters).

    Returns:
        Tuple[float, float]: (Easting, Northing) in UTM zone 27 coordinate system (meters).
    """
    N_f = X * 0.9996
    E_f = Y * 0.9996
    E_UTM27 = E_f + 500000.0
    N_UTM27 = N_f + 10000000.0
    print(f"E:{E_UTM27:.3f}, N:{N_UTM27:.3f}")
