import sys


# ---------- Core: Myers linear-space LCS on int lists ----------
# returns (pi, pj): matched index pairs, increasing order
def myers_pairs(a, b):
    pi = []
    pj = []

    def middle(a0, a1, b0, b1):
        N = a1 - a0
        M = b1 - b0
        sa = a[a0:a1]
        sb = b[b0:b1]
        ra = sa[::-1]
        rb = sb[::-1]
        # sentinels: stop snakes at the end without bounds checks
        sa.append(-1)
        ra.append(-1)
        sb.append(-2)
        rb.append(-2)
        delta = N - M
        odd = delta & 1
        max_d = (N + M + 1) // 2
        o = max_d + 1
        size = 2 * max_d + 3
        vf = [0] * size
        vb = [0] * size
        # active diagonal windows (shrink when a path runs off the grid)
        fs, fe = 0, 0
        bs, be = 0, 0
        for d in range(max_d + 1):
            # ---- forward ----
            for k in range(-d + fs, d - fe + 1, 2):
                if k == -d or (k != d and vf[o + k - 1] < vf[o + k + 1]):
                    x = vf[o + k + 1]
                else:
                    x = vf[o + k - 1] + 1
                x0 = x
                y = x - k
                if x < N and y < M:
                    while sa[x] == sb[y]:
                        x += 1
                        y += 1
                vf[o + k] = x
                if x > N:
                    fe += 2
                elif y > M:
                    fs += 2
                elif odd:
                    kk = delta - k
                    if -(d - 1) <= kk <= d - 1 and x + vb[o + kk] >= N:
                        return a0 + x0, b0 + x0 - k, a0 + x, b0 + y
            # ---- backward ----
            for k in range(-d + bs, d - be + 1, 2):
                if k == -d or (k != d and vb[o + k - 1] < vb[o + k + 1]):
                    x = vb[o + k + 1]
                else:
                    x = vb[o + k - 1] + 1
                x0 = x
                y = x - k
                if x < N and y < M:
                    while ra[x] == rb[y]:
                        x += 1
                        y += 1
                vb[o + k] = x
                if x > N:
                    be += 2
                elif y > M:
                    bs += 2
                elif not odd:
                    kk = delta - k
                    if -d <= kk <= d and x + vf[o + kk] >= N:
                        return (a0 + N - x, b0 + M - y,
                                a0 + N - x0, b0 + M - (x0 - k))
        raise RuntimeError("middle snake not found")

    def solve(a0, a1, b0, b1):
        while a0 < a1 and b0 < b1 and a[a0] == b[b0]:
            pi.append(a0)
            pj.append(b0)
            a0 += 1
            b0 += 1
        s = 0
        while a0 < a1 and b0 < b1 and a[a1 - 1] == b[b1 - 1]:
            a1 -= 1
            b1 -= 1
            s += 1
        if a0 < a1 and b0 < b1:
            sx, sy, ex, ey = middle(a0, a1, b0, b1)
            solve(a0, sx, b0, sy)
            for t in range(ex - sx):
                pi.append(sx + t)
                pj.append(sy + t)
            solve(ex, a1, ey, b1)
        for t in range(s):
            pi.append(a1 + t)
            pj.append(b1 + t)

    solve(0, len(a), 0, len(b))
    return pi, pj


# ---------- LCS with prefix/suffix trim + unique-line filter ----------
def lcs_pairs(ia, ib):
    n, m = len(ia), len(ib)
    lim = min(n, m)
    p = 0
    while p < lim and ia[p] == ib[p]:
        p += 1
    s = 0
    while s < lim - p and ia[n - 1 - s] == ib[m - 1 - s]:
        s += 1
    mi = list(range(p))
    mj = list(range(p))
    ma = ia[p:n - s]
    mb = ib[p:m - s]
    if ma and mb:
        in_a = set(ma)
        in_b = set(mb)
        ka = [i for i, v in enumerate(ma) if v in in_b]
        kb = [j for j, v in enumerate(mb) if v in in_a]
        ra = [ma[i] for i in ka]
        rb = [mb[j] for j in kb]
        pi, pj = myers_pairs(ra, rb)
        mi.extend([p + ka[x] for x in pi])
        mj.extend([p + kb[y] for y in pj])
    mi.extend(range(n - s, n))
    mj.extend(range(m - s, m))
    return mi, mj


def read_lines(path):
    with open(path, "rb") as f:
        data = f.read()
    parts = data.split(b"\n")
    if parts and parts[-1] == b"":
        parts.pop()
    return parts


def gaps(matched, length):
    r = []
    prev = 0
    for idx in matched:
        if idx > prev:
            r.append((prev, idx))
        prev = idx + 1
    if length > prev:
        r.append((prev, length))
    return r


def fmt(r):
    if not r:
        return "."
    return ",".join("%d-%d" % (s, e) for s, e in r)


def char_ranges(old, new):
    if old == new:
        return ".", "."
    mi, mj = lcs_pairs([ord(c) for c in old], [ord(c) for c in new])
    return fmt(gaps(mi, len(old))), fmt(gaps(mj, len(new)))


def render(A, B, highlight):
    table = {}
    ia = [table.setdefault(l, len(table)) for l in A]
    ib = [table.setdefault(l, len(table)) for l in B]
    mi, mj = lcs_pairs(ia, ib)
    mi.append(len(A))      # sentinel keep at the end
    mj.append(len(B))
    out = []
    pa = pb = 0
    for i, j in zip(mi, mj):
        nd = i - pa
        ni = j - pb
        if nd or ni:
            for t in range(pa, i):
                out.append(b"-" + A[t] + b"\n")
            for t in range(ni):
                out.append(b"+" + B[pb + t] + b"\n")
                if highlight and t < nd:
                    old = A[pa + t].decode("utf-8", "surrogateescape")
                    new = B[pb + t].decode("utf-8", "surrogateescape")
                    ro, rn = char_ranges(old, new)
                    out.append(("? %s | %s\n" % (ro, rn)).encode())
        if i < len(A):
            out.append(b" " + A[i] + b"\n")
        pa = i + 1
        pb = j + 1
    sys.stdout.buffer.write(b"".join(out))
    sys.stdout.buffer.flush()


def main():
    if len(sys.argv) != 4 or sys.argv[1] not in ("lines", "highlight"):
        sys.stderr.write("usage: main.py lines|highlight A B\n")
        sys.exit(2)
    try:
        A = read_lines(sys.argv[2])
        B = read_lines(sys.argv[3])
    except OSError as e:
        sys.stderr.write("error: %s\n" % e)
        sys.exit(2)
    render(A, B, sys.argv[1] == "highlight")


if __name__ == "__main__":
    main()