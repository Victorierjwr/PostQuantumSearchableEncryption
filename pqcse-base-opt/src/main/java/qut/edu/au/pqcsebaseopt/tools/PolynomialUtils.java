package qut.edu.au.pqcsebaseopt.tools;

import java.math.BigInteger;
import java.util.Arrays;

public class PolynomialUtils {
    /**
     * 根据根计算多项式系数 (Polynomial Expansion from Roots)
     * 计算 P(x) = (x - roots[0]) * ... * (x - roots[N-1])
     * * @param roots 根的数组 x_1, ..., x_N
     * @param q 模数
     * @return 系数数组 {a_0, a_1, ..., a_N}，其中 a_i 是 x^i 的系数
     */
    public static BigInteger[] getCoeffsFromRoots(BigInteger[] roots, BigInteger q) {
        int N = roots.length;
        BigInteger[] coeffs = new BigInteger[N + 1];

        // 初始化，使用 ZERO 填充，避免空指针
        Arrays.fill(coeffs, BigInteger.ZERO);
        coeffs[0] = BigInteger.ONE;

        for (int i = 0; i < N; i++) {
            BigInteger r = roots[i];

            // 1. 最高位特殊处理 (利用稀疏性)
            coeffs[i + 1] = coeffs[i];

            // 2. 中间位更新
            for (int j = i; j >= 1; j--) {
                // 原逻辑：c[j] = c[j-1] - r * c[j]
                coeffs[j] = coeffs[j - 1].subtract(coeffs[j].multiply(r)).mod(q);
            }

            // 3. 常数项更新
            // c[0] = - c[0] * r
            coeffs[0] = coeffs[0].multiply(r).negate().mod(q);
        }

        return coeffs;
    }
}
