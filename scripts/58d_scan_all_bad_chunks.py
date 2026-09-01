# 58d_scan_all_bad_chunks.py — 全量扫描 gzip/bgzf 文件的所有损坏区段, 输出坏块索引集合
# 原理: 顺序解析 member (BGZF member < 66KB, 窗口式喂数据避免大拷贝);
#       member 失败时向后搜 magic 重同步, 逐段报告损坏区间
# 用法: python 58d_scan_all_bad_chunks.py <file> <chunk_size>
import sys
import zlib

path = sys.argv[1]
CS = int(sys.argv[2])

BGZF_MAGIC = b"\x1f\x8b\x08\x04"
GZ_MAGIC = b"\x1f\x8b\x08"
WIN = 70000  # BGZF member 上限 ~66KB

f = open(path, "rb")
n = f.seek(0, 2)
f.seek(0)
print(f"file size: {n}")

bad_regions = []  # (start, next_good_start)
pos = 0

def try_member(pos):
    """在 pos 处尝试解压一个 member; 返回 ('ok', next_pos) / ('bad', None) / ('eofstream', None)"""
    f.seek(pos)
    window = f.read(WIN)
    if not window:
        return ("eofstream", None)
    try:
        d = zlib.decompressobj(zlib.MAX_WBITS | 16)
        d.decompress(window)
        if d.eof:
            # 注意: 窗口可能短于 WIN (文件尾), 按实际窗口长度计算
            return ("ok", pos + len(window) - len(d.unused_data))
        # 未到 member 尾: 若已到文件尾则截断, 否则窗口太小? (对 bgzf 不会发生)
        if pos + len(window) >= n:
            return ("bad", None)   # 截断的 member
        return ("bad", None)
    except zlib.error:
        return ("bad", None)

def resync(from_pos):
    """从 from_pos 向后找下一个能成功解压的 magic 起点"""
    f.seek(from_pos)
    tail = f.read(n - from_pos)
    q = 0
    while True:
        k = tail.find(BGZF_MAGIC, q)
        g = tail.find(GZ_MAGIC, q)
        cands = [x for x in (k, g) if x != -1]
        if not cands:
            return None
        c = min(cands)
        st, _ = try_member(from_pos + c)
        if st == "ok":
            return from_pos + c
        q = c + 1

while pos < n:
    st, nxt = try_member(pos)
    if st == "ok":
        pos = nxt
        continue
    if st == "eofstream":
        break
    # 损坏 member: 重同步
    found = resync(pos + 1)
    bad_regions.append((pos, found if found is not None else n))
    if found is None:
        break
    pos = found

if not bad_regions:
    print("ALL OK")
    sys.exit(0)
bad_chunks = set()
for s, e in bad_regions:
    print(f"corrupt region: [{s}, {e})")
    for ci in range(s // CS, min(e // CS + 1, (n + CS - 1) // CS)):
        bad_chunks.add(ci)
print("BAD_CHUNKS " + ",".join(str(i) for i in sorted(bad_chunks)))
