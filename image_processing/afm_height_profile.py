"""
AFM TIFF image analyzer: plots height vs position along the middle horizontal line.

Supports Park Systems (XEI/SmartScan) TIFF via custom tags 50434/50435,
Gwyddion key=value ImageDescription, Bruker NanoScope \\Key: Value format,
and standard TIFF XResolution/YResolution tags.

Falls back to prompting for any value not found in metadata.
Output is saved to image_processing/output/<stem>_profile.png.

Usage:
    python afm_height_profile.py
"""

import os
import re
import struct
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")

# Standard TIFF tag IDs
_TAG_IMAGE_DESC = 270
_TAG_X_RES      = 282
_TAG_Y_RES      = 283
_TAG_RES_UNIT   = 296

# Park Systems custom TIFF tag IDs
_TAG_PARK_DATA   = 50434   # PSIAData  – raw int16 height pixels
_TAG_PARK_HEADER = 50435   # PSIAHeader – binary parameter block

# Park PSIAHeader byte offsets (0-based; source: MATLAB AFMimage.m)
# Formats use struct notation: '<d' = little-endian double, '<i' = little-endian int32
_PK_NWIDTH        = (100, '<i')   # image width  (int32)
_PK_NHEIGHT       = (104, '<i')   # image height (int32)
_PK_X_SCAN_UM     = (140, '<d')   # X scan size  (double, µm)
_PK_Y_SCAN_UM     = (148, '<d')   # Y scan size  (double, µm)
_PK_DATA_GAIN     = (220, '<d')   # ADC-to-unit gain (double)
_PK_Z_SCALE       = (228, '<d')   # always 1.0
_PK_Z_OFFSET      = (236, '<d')   # always 0.0
_PK_UNIT_W        = (244, 16)     # UTF-16LE unit string (16 bytes → 8 chars)
_PK_DATA_MIN      = (260, '<i')   # min raw value (int32)
_PK_DATA_MAX      = (264, '<i')   # max raw value (int32)

_SI_FACTORS = {
    "pm": 1e-12, "nm": 1e-9, "µm": 1e-6, "um": 1e-6,
    "mm": 1e-3,  "cm": 1e-2, "m":  1.0,
}


def _to_meters(value: float, unit: str) -> float | None:
    f = _SI_FACTORS.get(unit.strip().lower().replace("\xb5", "µ"))
    return value * f if f is not None else None


def _parse_float(s) -> float | None:
    try:
        return float(str(s).strip())
    except (ValueError, TypeError):
        return None


def _rational(r) -> float | None:
    """Convert a Pillow TIFF rational to float."""
    if isinstance(r, (int, float)):
        return float(r)
    if isinstance(r, tuple):
        if len(r) == 2 and isinstance(r[0], int):
            return r[0] / r[1] if r[1] else None
        try:
            return r[0][0] / r[0][1]
        except (IndexError, TypeError, ZeroDivisionError):
            return None
    return None


def _unpack(header: bytes, offset: int, fmt):
    """Unpack a single value from header bytes; fmt is a struct format string."""
    size = struct.calcsize(fmt)
    if offset + size > len(header):
        return None
    return struct.unpack_from(fmt, header, offset)[0]


# ---------------------------------------------------------------------------
# Park Systems parser
# ---------------------------------------------------------------------------

def _parse_park_header(header: bytes) -> dict:
    """Decode the PSIAHeader binary block (tag 50435)."""
    if len(header) < 270:
        return {}

    result = {}

    width  = _unpack(header, *_PK_NWIDTH)
    height = _unpack(header, *_PK_NHEIGHT)
    if width and height and width > 0 and height > 0:
        result["park_width"]  = width
        result["park_height"] = height

    x_um = _unpack(header, *_PK_X_SCAN_UM)
    y_um = _unpack(header, *_PK_Y_SCAN_UM)
    if x_um and x_um > 0:
        result["x_range_m"] = x_um * 1e-6
    if y_um and y_um > 0:
        result["y_range_m"] = y_um * 1e-6

    gain   = _unpack(header, *_PK_DATA_GAIN)
    zscale = _unpack(header, *_PK_Z_SCALE) or 1.0
    zoff   = _unpack(header, *_PK_Z_OFFSET) or 0.0
    result["park_gain"]   = gain
    result["park_zscale"] = zscale
    result["park_zoff"]   = zoff

    data_min = _unpack(header, *_PK_DATA_MIN)
    data_max = _unpack(header, *_PK_DATA_MAX)
    result["park_data_min"] = data_min
    result["park_data_max"] = data_max

    unit_bytes = header[_PK_UNIT_W[0]: _PK_UNIT_W[0] + _PK_UNIT_W[1]]
    try:
        unit_str = unit_bytes.decode("utf-16-le").rstrip("\x00").strip()
    except UnicodeDecodeError:
        unit_str = unit_bytes.decode("ascii", errors="ignore").rstrip("\x00").strip()
    result["park_unit"] = unit_str if unit_str else "nm"

    return result


def _load_park_data(pil_img: Image.Image, meta: dict) -> np.ndarray | None:
    """
    Load calibrated height data from tag 50434.

    The standard TIFF image (pil_img.size) is often just a small embedded
    preview thumbnail with unrelated dimensions — the real scan resolution
    is given by the PSIAHeader's nWidth/nHeight fields (tag 50435). Sample
    width also varies by header version: some export float32 samples
    already normalised to the gain's units, others raw int16 counts. We
    pick whichever dtype's expected byte count exactly matches the tag's
    actual size, rather than guessing and silently truncating.

    Returns a 2-D float64 array in the units given by meta['park_unit'],
    or None if the tag is absent or its size doesn't match either dtype.
    """
    tags  = getattr(pil_img, "tag_v2", {})
    raw   = tags.get(_TAG_PARK_DATA)
    if raw is None:
        return None

    raw_bytes = bytes(raw) if not isinstance(raw, (bytes, bytearray)) else raw
    n_cols = meta.get("park_width")  or pil_img.size[0]
    n_rows = meta.get("park_height") or pil_img.size[1]
    n_pixels = n_cols * n_rows

    gain   = meta.get("park_gain") or 1.0
    zscale = meta.get("park_zscale") or 1.0
    zoff   = meta.get("park_zoff")   or 0.0

    if len(raw_bytes) == n_pixels * 4:
        z_raw = np.frombuffer(raw_bytes, dtype="<f4").reshape(n_rows, n_cols).astype(np.float64)
    elif len(raw_bytes) == n_pixels * 2:
        z_raw = np.frombuffer(raw_bytes, dtype="<i2").reshape(n_rows, n_cols).astype(np.float64)
    else:
        return None

    return gain * (zscale * z_raw + zoff)


# ---------------------------------------------------------------------------
# Gwyddion parser
# ---------------------------------------------------------------------------

def _parse_gwyddion(desc: str) -> dict:
    """Parse Gwyddion GWY TIFF ImageDescription (key=value, SI base units)."""
    kv = {}
    for line in desc.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            kv[k.strip().lower()] = v.strip()

    result = {}
    for meta_key, kv_key, unit_key in [
        ("x_range_m", "xreal", "xunit"),
        ("y_range_m", "yreal", "yunit"),
        ("z_range_m", "zreal", "zunit"),
    ]:
        val = _parse_float(kv.get(kv_key))
        if val is not None:
            unit  = kv.get(unit_key, "m")
            meters = _to_meters(val, unit)
            result[meta_key] = meters if meters is not None else val
    return result


# ---------------------------------------------------------------------------
# Bruker NanoScope parser
# ---------------------------------------------------------------------------

def _parse_nanoscope(desc: str) -> dict:
    """Parse Bruker/Veeco NanoScope \\Key: Value style metadata."""
    result = {}
    unit_pat = r"([\d.eE+\-]+)\s*(pm|nm|µm|um|mm|cm|m)\b"

    m = re.search(r"[Ss]can\s+[Ss]ize[:\s]+" + unit_pat, desc)
    if m:
        meters = _to_meters(float(m.group(1)), m.group(2))
        if meters:
            result["x_range_m"] = meters
            result["y_range_m"] = meters

    for pat in [r"[Zz]\s*[Rr]ange[:\s]+" + unit_pat,
                r"[Zz]\s*[Ss]cale[:\s]+" + unit_pat,
                r"[Hh]eight\s*[Rr]ange[:\s]+" + unit_pat]:
        m = re.search(pat, desc)
        if m:
            meters = _to_meters(float(m.group(1)), m.group(2))
            if meters:
                result["z_range_m"] = meters
                break
    return result


# ---------------------------------------------------------------------------
# Standard resolution tags
# ---------------------------------------------------------------------------

def _parse_resolution_tags(img: Image.Image, n_cols: int, n_rows: int) -> dict:
    """Use XResolution/YResolution + ResolutionUnit TIFF tags."""
    tags     = getattr(img, "tag_v2", {})
    res_unit = tags.get(_TAG_RES_UNIT, 2)
    unit_m   = {2: 0.0254, 3: 0.01}.get(res_unit)
    if unit_m is None:
        return {}

    result = {}
    for meta_key, tag_id, n_pix in [
        ("x_range_m", _TAG_X_RES, n_cols),
        ("y_range_m", _TAG_Y_RES, n_rows),
    ]:
        raw = tags.get(tag_id)
        if raw is not None:
            ppu = _rational(raw)
            if ppu and ppu > 0:
                result[meta_key] = (n_pix / ppu) * unit_m
    return result


# ---------------------------------------------------------------------------
# Top-level metadata reader
# ---------------------------------------------------------------------------

def read_afm_metadata(img: Image.Image, n_cols: int, n_rows: int) -> dict:
    """
    Try all known formats. Returns dict with keys:
      x_range_m, y_range_m, z_range_m  (metres, or None)
      park_gain, park_zscale, park_zoff, park_unit  (Park-specific, may be absent)
    """
    meta = {"x_range_m": None, "y_range_m": None, "z_range_m": None}

    tags = getattr(img, "tag_v2", {})

    # --- Park Systems ---
    park_header = tags.get(_TAG_PARK_HEADER)
    if park_header is not None:
        park_bytes = bytes(park_header) if not isinstance(park_header, (bytes, bytearray)) else park_header
        found = _parse_park_header(park_bytes)
        meta.update(found)
        # z_range derived from gain × data range; set a flag so caller knows
        meta["is_park"] = True
    else:
        meta["is_park"] = False

    # --- ImageDescription-based formats ---
    desc = tags.get(_TAG_IMAGE_DESC, "")
    if isinstance(desc, bytes):
        desc = desc.decode("utf-8", errors="ignore")

    for parser, args in [
        (_parse_gwyddion,       (desc,)),
        (_parse_nanoscope,      (desc,)),
        (_parse_resolution_tags, (img, n_cols, n_rows)),
    ]:
        found = parser(*args)
        for k, v in found.items():
            if meta.get(k) is None and v is not None:
                meta[k] = v

    return meta


# ---------------------------------------------------------------------------
# Image loading
# ---------------------------------------------------------------------------

def load_afm_tiff(path: str):
    """
    Returns (height_map, metadata, height_unit).
    Pixel data is loaded from the standard TIFF image channels (correct shape).
    For Park files the pixels are scaled by dfDataGain to give physical heights.
    """
    img    = Image.open(path)
    n_cols, n_rows = img.size
    meta   = read_afm_metadata(img, n_cols, n_rows)

    if meta.get("is_park") and meta.get("park_gain") is not None:
        park_data = _load_park_data(img, meta)
        if park_data is not None:
            return park_data, meta, meta.get("park_unit", "nm")

        # Fallback: tag 50434 (raw PSIAData) is missing, so reconstruct heights
        # by stretching the display TIFF's pixel range into gain * count range.
        # Less accurate than the raw counts above, since it assumes the
        # display image's grayscale span matches [data_min, data_max] exactly.
        data = np.array(img, dtype=np.float64)
        if data.ndim == 3:
            data = data[:, :, 0]

        gain     = meta["park_gain"]
        unit     = meta.get("park_unit", "nm")
        raw_min  = meta.get("park_data_min")
        raw_max  = meta.get("park_data_max")

        if raw_min is not None and raw_max is not None and raw_min != raw_max:
            z_range  = abs(gain) * abs(raw_max - raw_min)
            pix_min, pix_max = data.min(), data.max()
            data = z_range * (data - pix_min) / (pix_max - pix_min)
        else:
            data = -data * gain

        return data, meta, unit

    data = np.array(img, dtype=np.float64)
    if data.ndim == 3:
        data = data[:, :, 0]
    return data, meta, "raw counts"


# ---------------------------------------------------------------------------
# Profile extraction
# ---------------------------------------------------------------------------

def middle_horizontal_profile(height_map: np.ndarray):
    row = height_map.shape[0] // 2
    return row, height_map[row, :]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_profile(height_map: np.ndarray, profile: np.ndarray, row: int,
                 image_path: str, x_range_um: float | None,
                 z_unit: str, output_path: str):
    n_pixels = len(profile)

    baseline   = height_map.min()
    height_map = height_map - baseline
    profile    = profile - baseline

    if x_range_um is not None:
        x       = np.linspace(0, x_range_um, n_pixels)
        x_label = "Position (µm)"
    else:
        x       = np.arange(n_pixels)
        x_label = "Position (pixels)"

    z_label = f"Height ({z_unit})"

    # Slope: arctan(dz/dx) — normalise z and x to the same unit before atan
    dx = x[1] - x[0] if len(x) > 1 else 1.0
    dz_dx = np.gradient(profile, dx)
    if x_range_um is not None:
        z_per_um = _to_meters(1.0, z_unit)   # metres per z-unit
        if z_per_um is not None:
            dz_dx = dz_dx * (z_per_um / 1e-6)  # z-unit/µm → dimensionless
    slope_deg = np.degrees(np.arctan(dz_dx))

    fig, axes = plt.subplots(1, 3, figsize=(17, 4))

    ax_map = axes[0]
    im = ax_map.imshow(height_map, cmap="afmhot", aspect="equal",
                       interpolation="nearest")
    ax_map.axhline(row, color="cyan", linewidth=1.2, linestyle="--",
                   label=f"Row {row}")
    ax_map.set_title("AFM height map")
    ax_map.set_xlabel("x (pixels)")
    ax_map.set_ylabel("y (pixels)")
    ax_map.legend(fontsize=8)
    plt.colorbar(im, ax=ax_map, label=z_label)

    ax_prof = axes[1]
    ax_prof.plot(x, profile, linewidth=1.0, color="steelblue")
    ax_prof.set_title(f"Height profile — middle row ({row})")
    ax_prof.set_xlabel(x_label)
    ax_prof.set_ylabel(z_label)
    ax_prof.grid(True, linestyle=":", alpha=0.6)

    ax_slope = axes[2]
    ax_slope.plot(x, slope_deg, linewidth=1.0, color="darkorange")
    ax_slope.axhline(0, color="black", linewidth=0.6, linestyle="--")
    ax_slope.set_title(f"Slope — middle row ({row})")
    ax_slope.set_xlabel(x_label)
    ax_slope.set_ylabel("Slope (degrees)")
    ax_slope.grid(True, linestyle=":", alpha=0.6)

    fig.suptitle(os.path.basename(image_path), fontsize=9, y=1.01)
    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def ask(prompt: str) -> str:
    return input(f"{prompt}: ").strip()


def analyze_one():
    # Image path
    image_path = ""
    while not image_path:
        raw = ask("Path to AFM TIFF image").strip('"').strip("'")
        if os.path.isfile(raw):
            image_path = raw
        else:
            print(f"  File not found: {raw}")

    # Load
    height_map, meta, z_unit = load_afm_tiff(image_path)
    n_rows, n_cols = height_map.shape

    print(f"\nImage size : {n_cols} x {n_rows} px")
    if meta.get("is_park"):
        gain    = meta.get("park_gain", 0)
        dmin    = meta.get("park_data_min")
        dmax    = meta.get("park_data_max")
        z_range = abs(gain) * abs(dmax - dmin) if (dmin is not None and dmax is not None) else None
        print(f"Format     : Park Systems (PSIAHeader/PSIAData tags)")
        print(f"  gain = {gain:.6g} {z_unit}/count  |  raw counts: {dmin} to {dmax}")
        print(f"  z range  = {z_range:.4g} {z_unit}" if z_range is not None else "  z range  = unknown")
    else:
        print(f"Format     : generic TIFF / Gwyddion / NanoScope")

    def fmt_um(v_m):
        return f"{v_m*1e6:.4g} µm" if v_m is not None else "not found"
    def fmt_nm(v_m):
        return f"{v_m*1e9:.4g} nm" if v_m is not None else "not found"

    print(f"\nMetadata:")
    print(f"  X scan range : {fmt_um(meta['x_range_m'])}")
    print(f"  Y scan range : {fmt_um(meta['y_range_m'])}")
    if not meta.get("is_park"):
        print(f"  Z range      : {fmt_nm(meta['z_range_m'])}")

    # Prompt for missing X range
    if meta["x_range_m"] is None:
        s = ask("\nX scan range in µm (leave blank to use pixels)")
        meta["x_range_m"] = float(s) * 1e-6 if s else None

    x_range_um = meta["x_range_m"] * 1e6 if meta["x_range_m"] is not None else None

    row, profile = middle_horizontal_profile(height_map)
    print(f"\nSampled row : {row}")
    print(f"Height range on this row: {profile.min():.3f} – {profile.max():.3f} {z_unit}")

    stem        = os.path.splitext(os.path.basename(image_path))[0]
    output_path = os.path.join(OUTPUT_DIR, f"{stem}_profile.png")

    plot_profile(height_map, profile, row, image_path, x_range_um, z_unit, output_path)


def main():
    print("=== AFM Height Profile ===\n")

    while True:
        analyze_one()

        again = ask("\nScan another image? (y/n)").strip().lower()
        if again not in ("y", "yes"):
            break
        print()


if __name__ == "__main__":
    main()
