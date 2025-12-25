package qut.edu.au.pqcsebase.tools;

import java.math.BigInteger;

/**
 * 增量式线性无关检查器 (修正版)
 * 恢复了你的原始逻辑：如果 q=0，则使用默认大素数 1000000007 进行检查。
 */
public class IncrementalRankChecker {
    private final int dim;
    private final BigInteger q;
    private final BigInteger[][] basis;
    private final boolean[] pivotExists;

    // 你原来代码中的默认大素数
    private static final BigInteger DEFAULT_PRIME = BigInteger.valueOf(1000000007);

    public IncrementalRankChecker(int dim, BigInteger q) {
        this.dim = dim;
        // === 关键修复 ===
        // 如果传入 0 (BigInteger.ZERO)，则使用默认大素数
        // 这样 RandBasis 就可以利用它来做通用的线性无关检查，而不受格模数影响
        if (q == null || q.equals(BigInteger.ZERO)) {
            this.q = DEFAULT_PRIME;
        } else {
            this.q = q;
        }

        this.basis = new BigInteger[dim][dim];
        this.pivotExists = new boolean[dim];
    }

    public boolean checkAndAdd(BigInteger[] vec) {
        if (vec.length != dim) throw new IllegalArgumentException("Vector dimension mismatch");

        // 深拷贝并取模 (针对负数处理)
        BigInteger[] v = new BigInteger[dim];
        for (int i = 0; i < dim; i++) {
            v[i] = vec[i].mod(q);
            if (v[i].signum() < 0) v[i] = v[i].add(q);
        }

        for (int i = 0; i < dim; i++) {
            if (v[i].equals(BigInteger.ZERO)) continue;

            if (pivotExists[i]) {
                BigInteger pivot = basis[i][i];
                // 计算逆元 (在素数域下一定存在)
                BigInteger factor = v[i].multiply(pivot.modInverse(q)).mod(q);

                for (int j = i; j < dim; j++) {
                    BigInteger sub = factor.multiply(basis[i][j]).mod(q);
                    v[j] = v[j].subtract(sub).mod(q);
                    if (v[j].signum() < 0) v[j] = v[j].add(q);
                }
            } else {
                basis[i] = v;
                pivotExists[i] = true;
                return true;
            }
        }
        return false;
    }
}