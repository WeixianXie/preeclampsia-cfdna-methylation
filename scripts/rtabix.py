# -*- coding: utf-8 -*-
"""
rtabix.py - 纯 Python 远程 tabix 读取器 (HTTP Range + BGZF)
用于 GWAS Catalog FTP harmonised .h.tsv.gz (+.tbi) 的区域查询。
"""
import gzip
import struct
import subprocess
import time
import zlib


def http_get(url, start=None, end=None, timeout=120, retries=4):
    """GET 或 Range GET, 返回 bytes; 失败返回 None"""
    for i in range(retries):
        cmd = ["curl", "-s", "--ssl-no-revoke", "-m", str(timeout)]
        if start is not None:
            cmd += ["-r", "%d-%d" % (start, end)]
        cmd.append(url)
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=timeout + 15)
            if r.returncode == 0 and r.stdout:
                return r.stdout
        except Exception:
            pass
        time.sleep(2 + i * 2)
    return None


def _reg2bins(beg, end):
    """htslib reg2bins, 0-based [beg,end)"""
    lst = [0]
    end -= 1
    for k in range(1 + (beg >> 26), 2 + (end >> 26)):
        lst.append(k)
    for k in range(9 + (beg >> 23), 10 + (end >> 23)):
        lst.append(k)
    for k in range(73 + (beg >> 20), 74 + (end >> 20)):
        lst.append(k)
    for k in range(585 + (beg >> 17), 586 + (end >> 17)):
        lst.append(k)
    for k in range(4681 + (beg >> 14), 4682 + (end >> 14)):
        lst.append(k)
    return [b for b in lst if b < 4681 + 1 + (1 << 18)]


def _bgzf_block_size(buf, pos):
    """解析 BGZF 块 BC extra 字段, 返回块总字节数"""
    if pos + 18 > len(buf) or buf[pos] != 0x1F or buf[pos + 1] != 0x8B:
        return None
    xlen = struct.unpack_from("<H", buf, pos + 10)[0]
    p = pos + 12
    end_extra = pos + 12 + xlen
    while p + 4 <= min(end_extra, len(buf)):
        si1, si2, slen = buf[p], buf[p + 1], struct.unpack_from("<H", buf, p + 2)[0]
        if si1 == 66 and si2 == 67:  # 'BC'
            if p + 4 + 2 > len(buf):
                return None
            bsize = struct.unpack_from("<H", buf, p + 4)[0]
            return bsize + 1
        p += 4 + slen
    return None


class RemoteTabix:
    def __init__(self, url, cache_dir=None):
        self.url = url                      # .h.tsv.gz 完整 URL
        self.tbi_url = url + ".tbi"
        import os
        raw = None
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            cf = os.path.join(cache_dir, self.tbi_url.split("/")[-2] + "_" + self.tbi_url.split("/")[-1])
            if os.path.exists(cf) and os.path.getsize(cf) > 100:
                with open(cf, "rb") as fh:
                    raw = fh.read()
            else:
                raw = http_get(self.tbi_url)
                if raw:
                    with open(cf, "wb") as fh:
                        fh.write(raw)
        else:
            raw = http_get(self.tbi_url)
        if raw is None:
            raise IOError("无法下载 tbi: " + self.tbi_url)
        tbi = gzip.decompress(raw)
        self._parse_tbi(tbi)

    def _parse_tbi(self, t):
        magic, n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = \
            struct.unpack_from("<4siiiiiiii", t, 0)
        if magic != b"TBI\x01":
            raise ValueError("非 TBI 索引")
        self.n_ref = n_ref
        self.col_seq, self.col_beg, self.col_end = col_seq, col_beg, col_end
        self.meta, self.skip = meta, skip
        names_raw = t[36:36 + l_nm]
        self.refs = names_raw.decode().rstrip("\x00").split("\x00")
        p = 36 + l_nm
        self.bins = {}
        self.lins = {}
        for i in range(n_ref):
            (n_bin,) = struct.unpack_from("<i", t, p); p += 4
            binmap = {}
            for _ in range(n_bin):
                b, n_chunk = struct.unpack_from("<ii", t, p); p += 8
                chunks = []
                for _ in range(n_chunk):
                    cb, ce = struct.unpack_from("<QQ", t, p); p += 16
                    chunks.append((cb, ce))
                binmap[b] = chunks
            (n_intv,) = struct.unpack_from("<i", t, p); p += 4
            ioffs = struct.unpack_from("<%dQ" % n_intv, t, p); p += 8 * n_intv
            self.bins[i] = binmap
            self.lins[i] = ioffs
        # 尾部可能有 n_no_coor, 忽略

    def chrom_index(self, chrom):
        """chrom: 'chr6' 或 '6' -> 索引"""
        for i, nm in enumerate(self.refs):
            if nm == chrom or nm == chrom.replace("chr", "") or nm == "chr" + chrom:
                return i
        return None

    def _fetch_blocks(self, cbeg, cend):
        """取一个 chunk (虚拟偏移对) 的解压文本"""
        co0 = cbeg >> 16
        co1 = cend >> 16
        want_end = co1 + 65537
        buf = http_get(self.url, start=co0, end=want_end)
        if buf is None:
            return b""
        pos = 0
        out = bytearray()
        first = True
        while pos < len(buf):
            bsize = _bgzf_block_size(buf, pos)
            if bsize is None:
                break
            if pos + bsize > len(buf):
                # 需要更多数据
                more = http_get(self.url, start=co0 + len(buf), end=co0 + len(buf) + 65536)
                if not more:
                    break
                buf = buf + more
                continue
            d = zlib.decompressobj(31)
            try:
                block = d.decompress(bytes(buf[pos:pos + bsize]))
            except Exception:
                break
            cur_coff = co0 + pos
            if first:
                block = block[cbeg & 0xFFFF:]
                first = False
            if cur_coff == (cend >> 16):
                block = block[:cend & 0xFFFF]
                out += block
                break
            out += block
            pos += bsize
        return bytes(out)

    def fetch(self, chrom, beg, end):
        """0-based [beg,end) 区域查询, 返回行列表 (str)"""
        ci = self.chrom_index(chrom)
        if ci is None:
            return []
        b0, b1 = beg, end
        min_off = None
        for w in range(b0 >> 14, ((b1 - 1) >> 14) + 1):
            if w < len(self.lins[ci]) and self.lins[ci][w]:
                v = self.lins[ci][w]
                min_off = v if min_off is None else min(min_off, v)
        chunks = []
        for b in _reg2bins(b0, b1):
            for cb, ce in self.bins[ci].get(b, []):
                if min_off is not None and ce <= min_off:
                    continue
                chunks.append((cb, ce))
        if not chunks:
            return []
        chunks.sort()
        merged = []
        for cb, ce in chunks:
            if merged and cb <= merged[-1][1]:
                if ce > merged[-1][1]:
                    merged[-1] = (merged[-1][0], ce)
            else:
                merged.append((cb, ce))
        text = b""
        for cb, ce in merged:
            text += self._fetch_blocks(cb, ce)
        lines = []
        for ln in text.split(b"\n"):
            if not ln:
                continue
            s = ln.decode("utf-8", "replace")
            if s.startswith("#") or s.startswith("chr\ta"):
                continue
            lines.append(s)
        return lines

    def header(self):
        """取文件头 (列名行)"""
        buf = http_get(self.url, start=0, end=65535)
        if buf is None:
            return None
        pos = 0
        d = zlib.decompressobj(31)
        block = d.decompress(bytes(buf))
        for ln in block.split(b"\n"):
            s = ln.decode("utf-8", "replace")
            if s and not s.startswith("#"):
                return s
            if s.startswith("#") and len(s) > 2 and s[1] != "#":
                return s.lstrip("#")
        return None
