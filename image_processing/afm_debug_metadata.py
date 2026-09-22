"""
Diagnostic: dump raw Park AFM metadata + raw/display pixel stats for one file.

Usage:
    python afm_debug_metadata.py path/to/scan.tif
"""
import sys
import numpy as np
from PIL import Image

from afm_height_profile import (
    read_afm_metadata, _load_park_data, _TAG_PARK_DATA, _TAG_PARK_HEADER,
    _unpack, _PK_NWIDTH, _PK_NHEIGHT,
)


def main():
    path = sys.argv[1]
    img = Image.open(path)
    n_cols, n_rows = img.size
    meta = read_afm_metadata(img, n_cols, n_rows)

    tags = getattr(img, "tag_v2", {})
    header = tags.get(_TAG_PARK_HEADER)
    raw_tag = tags.get(_TAG_PARK_DATA)

    print(f"Image size (PIL)  : {n_cols} x {n_rows}")
    print(f"Raw tag 256/257   : ImageWidth={tags.get(256)} ImageLength={tags.get(257)}")
    print(f"RowsPerStrip(278) : {tags.get(278)}   StripByteCounts(279): {tags.get(279)}")
    print(f"BitsPerSample(258): {tags.get(258)}   SamplesPerPixel(277): {tags.get(277)}")
    print(f"PIL mode          : {img.mode}")
    print(f"Has 50435 header  : {header is not None} (len={len(header) if header else 0})")
    print(f"Has 50434 data    : {raw_tag is not None} (len={len(bytes(raw_tag)) if raw_tag else 0}, "
          f"expected={n_cols*n_rows*2})")
    if header is not None:
        hb = bytes(header)
        hw = _unpack(hb, *_PK_NWIDTH)
        hh = _unpack(hb, *_PK_NHEIGHT)
        print(f"Header nWidth/nHeight (offsets 100/104): {hw} x {hh}"
              + (f"  -> {hw*hh*2} bytes expected for 50434" if hw and hh else ""))
    print()
    print("Parsed metadata:")
    for k, v in meta.items():
        print(f"  {k:16s} = {v}")

    print()
    disp = np.array(img, dtype=np.float64)
    if disp.ndim == 3:
        disp = disp[:, :, 0]
    print(f"Display pixel data: min={disp.min():.6g} max={disp.max():.6g} dtype={img.mode}")

    if raw_tag is not None:
        z = _load_park_data(img, meta)
        if z is not None:
            print(f"Raw counts (int16): min={np.frombuffer(bytes(raw_tag)[:n_cols*n_rows*2], dtype='<i2').min()} "
                  f"max={np.frombuffer(bytes(raw_tag)[:n_cols*n_rows*2], dtype='<i2').max()}")
            print(f"Calibrated height  ({meta.get('park_unit')}): min={z.min():.6g} max={z.max():.6g} "
                  f"range={z.max()-z.min():.6g}")

    gain = meta.get("park_gain")
    dmin = meta.get("park_data_min")
    dmax = meta.get("park_data_max")
    if gain is not None and dmin is not None and dmax is not None:
        print(f"\nHeader gain={gain:.6g}, data_min={dmin}, data_max={dmax} "
              f"-> gain*(data_max-data_min) = {gain*(dmax-dmin):.6g} {meta.get('park_unit')}")

    if raw_tag is not None and header is not None:
        hb = bytes(header)
        hw = _unpack(hb, *_PK_NWIDTH)
        hh = _unpack(hb, *_PK_NHEIGHT)
        raw_bytes = bytes(raw_tag)
        if hw and hh:
            n = hw * hh
            print(f"\n--- int32 hypothesis: reshape as {hw} x {hh}, dtype=<i4 ---")
            i32 = np.frombuffer(raw_bytes[:n*4], dtype="<i4").reshape(hh, hw).astype(np.float64)
            print(f"raw int32: min={i32.min():.6g} max={i32.max():.6g}")
            z32 = gain * i32
            print(f"gain*raw (int32): min={z32.min():.6g} max={z32.max():.6g} range={z32.max()-z32.min():.6g} {meta.get('park_unit')}")

            print(f"\n--- float32 hypothesis: reshape as {hw} x {hh}, dtype=<f4 ---")
            f32 = np.frombuffer(raw_bytes[:n*4], dtype="<f4").reshape(hh, hw).astype(np.float64)
            print(f"raw float32 (as-is, no gain): min={f32.min():.6g} max={f32.max():.6g} range={f32.max()-f32.min():.6g}")
            print(f"gain*raw float32: min={(gain*f32).min():.6g} max={(gain*f32).max():.6g} range={(gain*f32).max()-(gain*f32).min():.6g}")


if __name__ == "__main__":
    main()
