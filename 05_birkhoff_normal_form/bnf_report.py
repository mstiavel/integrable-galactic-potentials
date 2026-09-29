import re, sys
for c in ['sph', 'flat', 'kk', 'core']:
    try: lines = open(f'bnfs_{c}.log').read().splitlines()
    except FileNotFoundError: continue
    hdr = [l for l in lines if l.startswith('---')]
    if not hdr: print(f"--- {c}: not finished"); continue
    print(hdr[0][:100])
    out = []
    for ln in lines:
        m = re.search(r'n=\s*(\d+)\s+max\|chi_n\|=([0-9.e+-]+)', ln)
        if m:
            n = int(m.group(1)); v = float(m.group(2))
            if n in (8, 10, 12, 14, 16, 18, 20): out.append(f"n={n}: {v**(1/n)/n:.4f}")
    print("   (max|chi_n|)^(1/n)/n : " + "  ".join(out))
