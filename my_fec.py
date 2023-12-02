#!/bin/python3

import typing

def bin_poly_mul(p1: int, p2: int) -> int:
    assert p1 >= 0
    assert p2 >= 0
    if p1 < p2:
        p1, p2 = p2, p1
    
    res = 0
    while p2 != 0:
        if (p2 & 1) != 0:
            res ^= p1
        p1 <<= 1
        p2 >>= 1
    
    return res

def bin_poly_div_mod(p: int, g: int) -> typing.Tuple[int, int]:
    assert p >= 0
    assert g > 0

    n = g.bit_length()
    k = p.bit_length()

    res = 0

    while k >= n:
        if (p & (1 << (k-1))) != 0:
            p ^= (g << (k-n))
            res |= (1<<(k-n))
        k -= 1
    
    return res, p

def bin_poly_mod(p: int, g: int) -> int:
    return bin_poly_div_mod(p, g)[1]

def bin_poly_mul_mod(p1: int, p2: int, g: int) -> int:
    return bin_poly_mod(bin_poly_mul(p1, p2), g)

def bin_poly_mul_mod(p1: int, p2: int, g: int) -> int:
    assert p1 >= 0
    assert p2 >= 0

    n = g.bit_length() - 1
    g -= (1<<n) # normalize g

    if p1 < p2:
        p1, p2 = p2, p1
    
    assert p1.bit_length() <= n

    res = 0
    while p2 != 0:
        if (p2 & 1) != 0:
            res ^= p1
        if (p1 & (1<<(n-1))) == 0:
            p1 <<= 1
        else:
            p1 = ((p1 << 1) & ((1<<n) - 1)) ^ g

        p2 >>= 1
    
    return res

def bin_poly_add(p1: int, p2: int) -> int:
    return p1 ^ p2

def bin_poly_euclid(p: int, g: int) -> typing.Tuple[typing.Tuple[int, int], int]:
    p_coefs, g_coefs = (1, 0), (0, 1)
    if p < g:
        p, g = g, p
        p_coefs, g_coefs = g_coefs, p_coefs
    
    # gcd(a,b)
    # a1 = b1*d1 + r1
    # a2 = b1, b2 = r1
    # a2 = b2*d2 + r2

    while g != 0:
        d, r = bin_poly_div_mod(p, g)
        # pi = gi*d + r
        # p_coef0*p + p_coef1*g = (g_coef0*p + g_coef1*g)*d + r
        # p*(p_coef0 - d*g_coef0) + g*(p_coef1 - d*g_coef1) = r
        # 
        p_coefs, g_coefs = g_coefs, (p_coefs[0] ^ bin_poly_mul(d, g_coefs[0]), p_coefs[1] ^ bin_poly_mul(d, g_coefs[1]))
        p, g = g, r

    return p_coefs, p

def bin_poly_inv_mod(p: int, g: int) -> int:
    return bin_poly_euclid(g, p)[0][1]

def bin_poly_pow_mod(p: int, n: int, g: int) -> int:
    res = 1
    cur_pow = p
    while n != 0:
        if (n & 1) != 0:
            res = bin_poly_mul_mod(res, cur_pow, g)
        cur_pow = bin_poly_mul_mod(cur_pow, cur_pow, g)
        n >>= 1
    return res

def bin_poly_inv_mod(p: int, g: int) -> int:
    n = g.bit_length() - 1
    return bin_poly_pow_mod(p, (1<<n) - 2, g)

def calc_cauchy_mat(n: int, k: int, g: int):
    assert n + k - 1 <= (1<<(g.bit_length() - 1))

    y = list(range(n))
    x = list(range(n, n+k-1))

    mat = [[1]*n]
    for i in range(k-1):
        row = []
        for j in range(n):
            row.append(bin_poly_inv_mod(bin_poly_add(x[i], y[j]), g))
        mat.append(row)
    
    return mat

def calc_inv_cauchy_mat(idxs: list[int], n: int, k:int, g: int):
    assert len(idxs) <= k
    has_one_row = True
    y = []
    x = list(range(n, n+k-1))
    for idx in idxs:
        if idx < n:
            y.append(idx)
        elif idx == n:
            has_one_row = False
        else:
            x.remove(idx - 1)
    
    if len(y) > 0:
        x = x[:len(y)-has_one_row]
    else:
        has_one_row = False
        x = []

    pi_xx = []
    for i in range(len(x)):
        res = 1
        for j in range(len(x)):
            if i == j:
                continue
            res = bin_poly_mul_mod(res, bin_poly_add(x[i], x[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_xx.append(res)
    
    pi_yy = []
    for i in range(len(y)):
        res = 1
        for j in range(len(y)):
            if i == j:
                continue
            res = bin_poly_mul_mod(res, bin_poly_add(y[i], y[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_yy.append(res)
    
    pi_xy = []
    for i in range(len(x)):
        res = 1
        for j in range(len(y)):
            res = bin_poly_mul_mod(res, bin_poly_add(x[i], y[j]), g)
        pi_xy.append(res)
    
    pi_yx = []
    for i in range(len(y)):
        res = 1
        for j in range(len(x)):
            res = bin_poly_mul_mod(res, bin_poly_add(y[i], x[j]), g)
        pi_yx.append(res)
    
    mat = []
    for i in range(len(y)):
        row = []
        for j in range(len(y)):
            if has_one_row:
                if j == 0:
                    row.append(bin_poly_mul_mod(pi_yx[i], pi_yy[i], g))
                    continue
                else:
                    j -= 1
            tmp = bin_poly_mul_mod(pi_xy[j], pi_yx[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_xx[j], g)
            tmp = bin_poly_mul_mod(tmp, pi_yy[i], g)
            tmp = bin_poly_mul_mod(tmp, bin_poly_inv_mod(bin_poly_add(x[j], y[i]), g), g)
            row.append(tmp)
        mat.append(row)
    return mat

def calc_inv_cauchy_mat2(idxs: list[int], n: int, k:int, g: int):
    assert len(idxs) <= k
    has_one_row = True
    y = []
    y_comp = list(range(n))
    x = list(range(n, n+k-1))
    for idx in idxs:
        if idx < n:
            y.append(idx)
            y_comp.remove(idx)
        elif idx == n:
            has_one_row = False
        else:
            x.remove(idx - 1)
    
    if len(y) > 0:
        x = x[:len(y)-has_one_row]
    else:
        has_one_row = False
        x = []
    
    assert len(y) + len(y_comp) == n

    pi_xx = []
    for i in range(len(x)):
        res = 1
        for j in range(len(x)):
            if i == j:
                continue
            res = bin_poly_mul_mod(res, bin_poly_add(x[i], x[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_xx.append(res)
    
    pi_yy = []
    for i in range(len(y)):
        res = 1
        for j in range(len(y)):
            if i == j:
                continue
            res = bin_poly_mul_mod(res, bin_poly_add(y[i], y[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_yy.append(res)
    
    pi_xy = []
    for i in range(len(x)):
        res = 1
        for j in range(len(y)):
            res = bin_poly_mul_mod(res, bin_poly_add(x[i], y[j]), g)
        pi_xy.append(res)
    
    pi_yx = []
    for i in range(len(y)):
        res = 1
        for j in range(len(x)):
            res = bin_poly_mul_mod(res, bin_poly_add(y[i], x[j]), g)
        pi_yx.append(res)
    
    pi_y_comp_x = []
    for i in range(len(y_comp)):
        res = 1
        for j in range(len(x)):
            res = bin_poly_mul_mod(res, bin_poly_add(y_comp[i], x[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_y_comp_x.append(res)
    
    pi_y_comp_y = []
    for i in range(len(y_comp)):
        res = 1
        for j in range(len(y)):
            res = bin_poly_mul_mod(res, bin_poly_add(y_comp[i], y[j]), g)
        pi_y_comp_y.append(res)
    
    mat = []
    for i in range(len(y)):
        row = []
        for j in range(len(y)):
            if has_one_row:
                if j == 0:
                    row.append(bin_poly_mul_mod(pi_yx[i], pi_yy[i], g))
                    continue
                else:
                    j -= 1
            tmp = bin_poly_mul_mod(pi_xy[j], pi_yx[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_xx[j], g)
            tmp = bin_poly_mul_mod(tmp, pi_yy[i], g)
            tmp = bin_poly_mul_mod(tmp, bin_poly_inv_mod(bin_poly_add(x[j], y[i]), g), g)
            row.append(tmp)
        mat.append(row)
    
    mat2 = []
    for i in range(len(y)):
        row = []
        for j in range(len(y_comp)):
            tmp = bin_poly_mul_mod(pi_y_comp_y[j], pi_yx[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_yy[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_y_comp_x[j], g)
            tmp = bin_poly_mul_mod(tmp, bin_poly_inv_mod(bin_poly_add(y_comp[j], y[i]), g), g)
            row.append(tmp)
        mat2.append(row)

    return mat, mat2



def calc_inv_cauchy_mat3(idxs: list[int], n: int, k:int, g: int):
    assert len(idxs) <= k
    has_one_row = True
    y = []
    y_comp = list(range(n))
    x = list(range(n, n+k-1))
    for idx in idxs:
        if idx < n:
            y.append(idx)
            y_comp.remove(idx)
        elif idx == n:
            has_one_row = False
        else:
            x.remove(idx - 1)
    
    if len(y) > 0:
        x = x[:len(y)-has_one_row]
    else:
        has_one_row = False
        x = []
    
    assert len(y) + len(y_comp) == n

    pi_xx = []
    for i in range(len(x)):
        res = 1
        for j in range(len(x)):
            if i == j:
                continue
            res = bin_poly_mul_mod(res, bin_poly_add(x[i], x[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_xx.append(res)
    
    pi_yy = []
    for i in range(len(y)):
        res = 1
        for j in range(len(y)):
            if i == j:
                continue
            res = bin_poly_mul_mod(res, bin_poly_add(y[i], y[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_yy.append(res)
    
    pi_xy = []
    for i in range(len(x)):
        res = 1
        for j in range(len(y)):
            res = bin_poly_mul_mod(res, bin_poly_add(x[i], y[j]), g)
        pi_xy.append(res)
    
    pi_yx = []
    for i in range(len(y)):
        res = 1
        for j in range(len(x)):
            res = bin_poly_mul_mod(res, bin_poly_add(y[i], x[j]), g)
        pi_yx.append(res)
    
    pi_y_comp_x = []
    for i in range(len(y_comp)):
        res = 1
        for j in range(len(x)):
            res = bin_poly_mul_mod(res, bin_poly_add(y_comp[i], x[j]), g)
        res = bin_poly_inv_mod(res, g)
        pi_y_comp_x.append(res)
    
    pi_y_comp_y = []
    for i in range(len(y_comp)):
        res = 1
        for j in range(len(y)):
            res = bin_poly_mul_mod(res, bin_poly_add(y_comp[i], y[j]), g)
        pi_y_comp_y.append(res)
    
    mat = []
    for i in range(len(y)):
        row = []
        for j in range(len(y)):
            if has_one_row:
                if j == 0:
                    row.append(bin_poly_mul_mod(pi_yx[i], pi_yy[i], g))
                    continue
                else:
                    j -= 1
            tmp = bin_poly_mul_mod(pi_xy[j], pi_yx[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_xx[j], g)
            tmp = bin_poly_mul_mod(tmp, pi_yy[i], g)
            tmp = bin_poly_mul_mod(tmp, bin_poly_inv_mod(bin_poly_add(x[j], y[i]), g), g)
            row.append(tmp)
        mat.append(row)
    
    mat2 = []
    for i in range(len(y)):
        row = []
        for j in range(len(y_comp)):
            tmp = bin_poly_mul_mod(pi_y_comp_y[j], pi_yx[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_yy[i], g)
            tmp = bin_poly_mul_mod(tmp, pi_y_comp_x[j], g)
            tmp = bin_poly_mul_mod(tmp, bin_poly_inv_mod(bin_poly_add(y_comp[j], y[i]), g), g)
            row.append(tmp)
        mat2.append(row)

    return mat, mat2


def mat_mul_vec(mat: list[list[int]], vec: list[int], g: int) -> list[int]:
    num_rows = len(mat)
    if num_rows == 0:
        return []
    num_cols = len(mat[0])
    assert len(vec) == num_cols

    res : list[int] = []
    for i in range(num_rows):
        res_i = 0
        for j in range(num_cols):
            res_i = bin_poly_add(res_i, bin_poly_mul_mod(mat[i][j], vec[j], g))
        res.append(res_i)
    return res

def mat_mul_mat(A: list[list[int]], B: list[list[int]], g: int):
    A_rows = len(A)
    if A_rows == 0:
        return []
    A_cols = len(A[0])
    B_rows = len(B)
    B_cols = len(B[0])
    if B_cols == 0:
        return []
    assert A_cols == B_rows

    res : list[list[int]] = []
    for i in range(A_rows):
        row = []
        for j in range(B_cols):
            res_ij = 0
            for k in range(A_cols):
                res_ij = bin_poly_add(res_ij, bin_poly_mul_mod(A[i][k], B[k][j], g))
            row.append(res_ij)
        res.append(row)
    return res

def my_coder_redundency(seq: list[int], k:int, g: int):
    n = len(seq)
    C = calc_cauchy_mat(n, k, g)
    return mat_mul_vec(C, seq, g)

def my_coder(seq: list[int], k:int, g: int):
    res = list(seq)
    res.extend(my_coder_redundency(seq, k, g))
    return res

def calc_cauchy_mat2(idxs: list[int], n: int, k:int, g: int):
    assert len(idxs) <= k
    has_one_row = True
    y = []
    y_comp = list(range(n))
    x = list(range(n, n+k-1))
    for idx in idxs:
        if idx < n:
            y.append(idx)
            y_comp.remove(idx)
        elif idx == n:
            has_one_row = False
        else:
            x.remove(idx - 1)
    
    if len(y) > 0:
        x = x[:len(y)-has_one_row]
    else:
        has_one_row = False
        x = []
    
    assert len(y) + len(y_comp) == n

    if has_one_row:
        mat = [[1]*len(y)]
    else:
        mat = []
    for i in range(len(x)):
        row = []
        for j in range(len(y)):
            row.append(bin_poly_inv_mod(bin_poly_add(x[i], y[j]), g))
        mat.append(row)
    
    if has_one_row:
        mat2 = [[1]*len(y_comp)]
    else:
        mat2 = []
    for i in range(len(x)):
        row = []
        for j in range(len(y_comp)):
            row.append(bin_poly_inv_mod(bin_poly_add(x[i], y_comp[j]), g))
        mat2.append(row)

    return mat, mat2

def my_decoder(seq: list[typing.Union[int, None]], k:int, g: int):
    assert len(seq) >= k
    n  = len(seq) - k
    present_seq = []
    idxs = []
    for i in range(len(seq)):
        x = seq[i]
        if x is not None:
            present_seq.append(x)
        else:
            idxs.append(i)
    
    assert len(idxs) <= k

    present_seq = present_seq[:n]
    C1_inv,C1_inv_C2 = calc_inv_cauchy_mat2(idxs, n, k, g)

    mat = [C1_inv_C2[i] + C1_inv[i] for i in range(len(C1_inv))]

    recovered = mat_mul_vec(mat, present_seq, g)

    res_seq : list[int] = []
    j = 0
    for i in range(n):
        x = seq[i]
        if x is not None:
            res_seq.append(x)
        else:
            res_seq.append(recovered[j])
            j += 1

    return res_seq

def main():
    g = (1<<16) + (1<<5) + (1<<3) + (1<<1) + (1<<0)
    p = (1<<0)
    
    C = calc_cauchy_mat(20, 10, g)
    C1_inv,C1_inv_C2 = calc_inv_cauchy_mat2([0,2,4], 20, 10, g)
    C1, C2 = calc_cauchy_mat2([0,2,4], 20, 10, g)
    print(C)
    print(C1)
    print(C2)
    print(C1_inv)
    print(C1_inv_C2)
    print(mat_mul_mat(C1, C1_inv, g))
    print(mat_mul_mat(C1, C1_inv_C2, g) == C2)


    seq = my_coder(list(range(10)), 3, g)
    print(seq)
    seq[2] = None
    seq[9] = None
    seq[6] = None
    print(seq)
    print(my_decoder(seq, 3, g))

    print(mat_mul_mat(C1, C1_inv_C2, g))
    

if __name__ == '__main__':
    main()
