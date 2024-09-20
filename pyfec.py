#!/usr/bin/env python3

import sys
import array
import typing

_is_sys_big_endian = (sys.byteorder == 'big')

POLY_G = (1<<5) + (1<<3) + (1<<1) + (1<<0)

def poly_add(a: int, b: int) -> int:
    return a ^ b

def poly_mul(a: int, b: int) -> int:
    res = 0
    for i in range(16):
        res <<= 1
        res ^= ((1<<16) | POLY_G) & (-(res >> 16))
        res ^= a & (-(b >> 15))
        b = (b << 1) & ((1<<16) - 1)
    
    return res

def poly_pow(a: int, n: int) -> int:
    res = 1
    a_2_i = a
    while n != 0:
        if (n & 1) != 0:
            res = poly_mul(res, a_2_i)
        a_2_i = poly_mul(a_2_i, a_2_i)
        n >>= 1
    
    return res

def poly_inv(a: int) -> int:
    return poly_pow(a, (1<<16) - 1 - 1)

class FecInvCache:
    def __init__(self, n_k_1: int) -> None:
        assert n_k_1 <= (1<<16)

        if n_k_1 > 1:
            sz = 1 << (n_k_1 - 1).bit_length()
        else:
            sz = 1
        
        self.inv_arr = [poly_inv(i) for i in range(1, sz)]
    
    @staticmethod
    def from_n_k(n: int, k: int) -> 'FecInvCache':
        assert n >= 0 and k >= 0
        # this is python so negative here is fine
        return FecInvCache(n + k - 1)
    
    def inv(self, a: int) -> int:
        assert a > 0
        return self.inv_arr[a - 1]

def fec_get_redundancy_pak(inv_cache: FecInvCache, paks: typing.Sequence[typing.Union[bytes, bytearray]], idx: int) -> bytes:
    i = 0
    n = len(paks)
    if n == 0:
        return b''
    
    assert idx >= 0
    assert n + idx <= (1<<16)
    assert idx == 0 or n + idx - 1 <= len(inv_cache.inv_arr)

    for pak in paks:
        if i == 0:
            l = len(pak)
            assert (l & 1) == 0
            res_pak = array.array('H', bytearray(l))
            l >>= 1
        else:
            assert (len(pak) == (l << 1))
        
        pak = array.array('H', pak)
        if _is_sys_big_endian:
            pak.byteswap()

        if idx == 0:
            for j in range(l):
                res_pak[j] ^= pak[j]
        else:
            a_i = inv_cache.inv(poly_add(n + idx - 1, i))
            for j in range(l):
                res_pak[j] ^= poly_mul(pak[j], a_i)
        
        i += 1
    
    if _is_sys_big_endian:
        res_pak.byteswap()

    return res_pak.tobytes()
    


class FecRxCtx:
    def __init__(self, n: int) -> None:
        assert n <= (1<<16) + 1
        assert n >= 0
        self.n = n
        self.l = None
        self.paks : typing.Dict[int, typing.Union[bytes, bytearray]] = dict()
    
    def add_pak(self, pak: typing.Union[bytes, bytearray], idx: int) -> bool:
        assert idx <= (1<<16)

        if self.l is None:
            assert (len(pak) & 1) == 0
            self.l = len(pak) >> 1
        else:
            assert len(pak) == (self.l << 1)
        
        if len(self.paks) < self.n and idx not in self.paks:
            self.paks[idx] = pak
        
        return self.can_reconstruct()
    
    def can_reconstruct(self) -> bool:
        return len(self.paks) == self.n

    def reconstruct(self, inv_cache: FecInvCache) -> bytes:
        assert self.can_reconstruct()

        if self.l is None: # which means n is 0
            return b''
        
        ones_pak = None
        x_paks : typing.List['array.array[int]'] = []
        y_paks : typing.List['array.array[int]'] = []
        n = self.n
        l = self.l

        x_arr : typing.List[int] = []
        y_arr : typing.List[int] = []
        max_x = n - 1

        for idx in self.paks:
            pak = array.array('H', self.paks[idx])
            if _is_sys_big_endian:
                pak.byteswap()

            if idx < n:
                y_paks.append(pak)
                y_arr.append(idx)
            elif idx == self.n:
                ones_pak = pak
            else:
                x_paks.append(pak)
                x_arr.append(idx - 1)
                max_x = max(max_x, idx - 1)
        
        assert len(inv_cache.inv_arr) >= max_x

        y_missing = [idx for idx in range(n) if idx not in self.paks]

        res_buf = array.array('H', bytearray(n*l*2))
        for i in range(len(y_paks)):
            y_i = y_arr[i]
            res_buf[y_i*l:(y_i+1)*l] = y_paks[i]
        
        if len(y_paks) == n:
            if _is_sys_big_endian:
                res_buf.byteswap()
            return res_buf.tobytes()

        for i in range(len(y_paks)):
            y_i = y_arr[i]
            pi_ycomp_y_div_ycomp_x_i = 1
            for j in range(len(y_missing)):
                pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, poly_add(y_i, y_missing[j]))
            for j in range(len(x_arr)):
                pi_ycomp_y_div_ycomp_x_i = poly_mul(pi_ycomp_y_div_ycomp_x_i, inv_cache.inv(poly_add(y_i, x_arr[j])))
            
            pak = y_paks[i]
            for ii in range(l):
                pak[ii] = poly_mul(pak[ii], pi_ycomp_y_div_ycomp_x_i)

        for i in range(len(x_paks)):
            x_i = x_arr[i]
            pi_xy_div_xx_i = 1
            for j in range(len(y_missing)):
                pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, poly_add(x_i, y_missing[j]))
            for j in range(len(x_arr)):
                if j != i:
                    pi_xy_div_xx_i = poly_mul(pi_xy_div_xx_i, inv_cache.inv(poly_add(x_i, x_arr[j])))
            
            pak = x_paks[i]
            for ii in range(l):
                pak[ii] = poly_mul(pak[ii], pi_xy_div_xx_i)
        
        for ii in range(l):
            ones_pak_ii = 0
            if ones_pak is not None:
                ones_pak_ii = ones_pak[ii]
            
            recovered = array.array('H', [ones_pak_ii]*len(y_missing))

            for j in range(len(y_paks)):
                y_j = y_arr[j]
                pak_val = y_paks[j][ii]
                for i in range(len(y_missing)):
                    y_i = y_missing[i]
                    recovered[i] ^= poly_mul(pak_val, inv_cache.inv(poly_add(y_i, y_j)))
            
            for j in range(len(x_paks)):
                x_j = x_arr[j]
                pak_val = x_paks[j][ii]
                for i in range(len(y_missing)):
                    y_i = y_missing[i]
                    recovered[i] ^= poly_mul(pak_val, inv_cache.inv(poly_add(y_i, x_j)))
            
            for i in range(len(y_missing)):
                res_buf[y_missing[i]*l + ii] = recovered[i]
        
        for i in range(len(y_missing)):
            y_i = y_missing[i]
            pi_yx_div_yy_i = 1
            for j in range(len(x_arr)):
                pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, poly_add(y_i, x_arr[j]))
            for j in range(len(y_missing)):
                if j != i:
                    pi_yx_div_yy_i = poly_mul(pi_yx_div_yy_i, inv_cache.inv(poly_add(y_i, y_missing[j])))
            
            for ii in range(l):
                res_buf[y_i*l + ii] = poly_mul(res_buf[y_i*l + ii], pi_yx_div_yy_i)
        
        if _is_sys_big_endian:
            res_buf.byteswap()
        return res_buf.tobytes()


def main():
    n = 100
    k = 10
    paks = [bytes(i for j in range(4)) for i in range(n)]
    inv_cache = FecInvCache.from_n_k(n, k)
    print(len(inv_cache.inv_arr))
    paks += [fec_get_redundancy_pak(inv_cache, paks, i) for i in range(k)]
    rx_ctx = FecRxCtx(n)
    for i in range(k, len(paks)):
        rx_ctx.add_pak(paks[i], i)
    print(repr(rx_ctx.reconstruct(inv_cache)))
    pass

if __name__ == '__main__':
    main()

