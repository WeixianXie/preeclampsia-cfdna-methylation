# 58c_find_bad_chunk.py — 定位 gzip/bgzip 文件中第一个损坏字节的压缩偏移, 映射到分块索引
# 用法: python 58c_find_bad_chunk.py <file> <chunk_size>
import sys
import zlib

path = sys.argv[1]
CS = int(sys.argv[2])

f = open(path, "rb")
d = zlib.decompressobj(zlib.MAX_WBITS | 16)
buf = b""
fpos = 0
bad = None

while True:
    chunk = f.read(1 << 20)
    if not chunk:
        break
    buf += chunk
    start_buf = fpos - len(buf) + len(buf)  # unused
    while buf:
        if d.eof:
            d = zlib.decompressobj(zlib.MAX_WBITS | 16)
        # 当前 member 从 (fpos - len(buf)) 开始
        member_start = fpos - len(buf)
        try:
            d.decompress(buf)
        except zlib.error as e:
            bad = member_start
            print(f"CORRUPT: member starting at compressed offset {member_start}: {e}")
            break
        if d.eof:
            buf = d.unused_data
        else:
            buf = b""
    if bad is not None:
        break
    fpos = f.tell()

if bad is None:
    print("ALL OK")
    sys.exit(0)
idx = bad // CS
print(f"损坏首个成员偏移 {bad} -> 分块索引 {idx} (part_{idx:05d}), 建议重下该分块及前后各一块")
