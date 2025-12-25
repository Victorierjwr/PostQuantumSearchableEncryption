package qut.edu.au.pqcsebase.algorithms;

import qut.edu.au.pqcsebase.tools.BigIntMatrix;
import qut.edu.au.pqcsebase.tools.IncrementalRankChecker;

import java.math.BigInteger;
import java.util.Arrays;

import static qut.edu.au.pqcsebase.algorithms.SampleD.*;
import static qut.edu.au.pqcsebase.tools.BigIntMatrix.*;

public class RandBasis {

    /**
     * Generate a random basis for the lattice defined by the basis S.
     *
     * @param S      The original basis matrix.
     * @param sigma      The Gaussian parameter for sampling.
     * @param nParam lattice A dimension nParam = n.
     * @return A new basis matrix sampled from the lattice defined by S.
     */
    public static BigIntMatrix randBasis(BigIntMatrix S, double sigma, long nParam) {
        int m = S.getRowDimension();
        int n = S.getColumnDimension();
        BigInteger q = S.getModulus();

        BigIntMatrix S_prime = new BigIntMatrix(m, n, q);
        double[] center = new double[m];
        Arrays.fill(center, 0.0);

        // 初始化优化后的检查器
        IncrementalRankChecker rankChecker = new IncrementalRankChecker(Math.max(n, m), BigInteger.ZERO);

        for (int i = 0; i < n; i++) {
            BigInteger[] sampledVector = sampleD(S, sigma, center, nParam);

            // 简单的零向量检查
            boolean isZero = true;
            for (BigInteger val : sampledVector) {
                if (!val.equals(BigInteger.ZERO)) {
                    isZero = false;
                    break;
                }
            }
            if (isZero) {
                i--; continue;
            }

            // 使用增量检查
            if (rankChecker.checkAndAdd(sampledVector)) {
                // 线性无关，填入矩阵
                for (int j = 0; j < m; j++) {
                    S_prime.set(j, i, sampledVector[j]);
                }
            } else {
                // 线性相关，重试
                i--;
            }
        }

        return LLL(S_prime, 0.99);
    }

    public static BigIntMatrix randBasisPre(BigIntMatrix S, double sigma, long nParam) {
        int m = S.getRowDimension();
        int n = S.getColumnDimension();
        BigInteger q = S.getModulus();

        BigIntMatrix S_prime = new BigIntMatrix(m, n, q);
        double[] center = new double[m];
        Arrays.fill(center, 0.0);

        double[][] B_double = new double[n][m];
        for(int i=0;i<n;i++) for(int j=0;j<m;j++) B_double[i][j] = S.get(j, i).doubleValue();
        double[][] B_tilde = new double[n][m];
        double[] normsSq = new double[n];
        computeGramSchmidt(B_double, B_tilde, normsSq);

        // 初始化优化后的检查器
        IncrementalRankChecker rankChecker = new IncrementalRankChecker(Math.max(n, m), BigInteger.ZERO);

        for (int i = 0; i < n; i++) {
            BigInteger[] sampledVector = sampleD_Fast(B_tilde, normsSq, S, sigma, center, nParam);

            boolean isZero = true;
            for (BigInteger val : sampledVector) {
                if (!val.equals(BigInteger.ZERO)) {
                    isZero = false;
                    break;
                }
            }
            if (isZero) {
                i--; continue;
            }

            // 使用增量检查
            if (rankChecker.checkAndAdd(sampledVector)) {
                for (int j = 0; j < m; j++) {
                    S_prime.set(j, i, sampledVector[j]);
                }
            } else {
                i--;
            }
            if(i % 100 == 0)System.out.println("Generated column " + (i + 1) + " / " + n);
        }

        return LLL(S_prime, 0.99);
    }

    /**
     * Lenstra–Lenstra–Lovász (LLL) lattice basis reduction algorithm.
     * Implements the requirements of "Lemma 1" by producing a delta-reduced basis.
     *
     * @param B     The input lattice basis (Column-Oriented).
     * @param delta The reduction parameter (typically 0.99 or 0.75).
     * @return A new BigIntMatrix containing the reduced basis.
     */
    public static BigIntMatrix LLL(BigIntMatrix B, double delta) {
        int n = B.getColumnDimension();
        int m = B.getRowDimension();

        // 1. Deep copy
        BigIntMatrix S = B.copy();

        // 2. GSO Cache
        // We only need mu and normsSq for the logic.
        // B_tilde vectors are NOT needed for the LLL logic itself if mu is maintained.
        double[] normsSq = new double[n];
        double[][] mu = new double[n][n];

        // Initial GSO (O(n^3)) - Done only once!
        computeGSO(S, normsSq, mu);

        int k = 1;
        while (k < n) {
            // === Step A: Size Reduction (Optimized) ===
            // Use cached mu table instead of computing dot products
            for (int j = k - 1; j >= 0; j--) {
                if (Math.abs(mu[k][j]) > 0.5) {
                    long c = Math.round(mu[k][j]);
                    if (c != 0) {
                        // Matrix operation: Col_k = Col_k - c * Col_j
                        S.colSub(k, j, BigInteger.valueOf(c));

                        // Update GSO coefficients (O(n) speed)
                        // mu[k][j] becomes small
                        mu[k][j] -= c;
                        // Update other coefficients affected by this subtraction
                        for (int l = 0; l < j; l++) {
                            mu[k][l] -= c * mu[j][l];
                        }
                    }
                }
            }

            // === Step B: Lovasz Condition ===
            double lhs = delta * normsSq[k - 1];
            double rhs = normsSq[k] + mu[k][k - 1] * mu[k][k - 1] * normsSq[k - 1];

            if (rhs < lhs) {
                // === Swap Columns k and k-1 ===
                S.swapCols(k, k - 1);

                // === Local GSO Update (The Magic Part) ===
                // Instead of recomputing the whole GSO (O(n^3)), we use algebra
                // to update only the changed values in O(n).

                int km1 = k - 1;
                double mu_k_km1 = mu[k][km1];
                double B_km1 = normsSq[km1]; // Old norm at k-1
                double B_k = normsSq[k];     // Old norm at k (orthogonal part)

                // 1. Update Norms
                double new_B_km1 = B_k + mu_k_km1 * mu_k_km1 * B_km1;
                double new_B_k = B_k * B_km1 / new_B_km1;

                normsSq[km1] = new_B_km1;
                normsSq[k] = new_B_k;

                // 2. Update mu[k][k-1]
                mu[k][km1] = mu_k_km1 * B_km1 / new_B_km1;

                // 3. Update mu matrix for columns j < k-1 (just swap)
                for (int j = 0; j < km1; j++) {
                    double tmp = mu[km1][j];
                    mu[km1][j] = mu[k][j];
                    mu[k][j] = tmp;
                }

                // 4. Update mu matrix for rows i > k (Givens rotation-like update)
                for (int i = k + 1; i < n; i++) {
                    double t = mu[i][k];
                    double u = mu[i][km1];
                    // These formulas preserve the validity of the GSO without recomputing vectors
                    mu[i][k] = u - mu_k_km1 * t;
                    mu[i][km1] = t + mu[k][km1] * mu[i][k];
                }

                // Backtrack
                k = Math.max(k - 1, 1);
            } else {
                k++;
            }
        }
        return S;
    }

    /**
     * Initial GSO computation.
     * Computes mu and normsSq from scratch.
     */
    private static void computeGSO(BigIntMatrix S, double[] normsSq, double[][] mu) {
        int n = S.getColumnDimension();
        int m = S.getRowDimension();

        // Temporary storage for orthogonal vectors (only needed during init)
        double[][] b_tilde = new double[n][m];

        for (int i = 0; i < n; i++) {
            // Load B[i] into b_tilde[i]
            for (int r = 0; r < m; r++) {
                b_tilde[i][r] = S.get(r, i).doubleValue();
            }

            // Orthogonalize against previous
            for (int j = 0; j < i; j++) {
                double dot = 0;
                for (int r = 0; r < m; r++) dot += b_tilde[i][r] * b_tilde[j][r];

                mu[i][j] = dot / normsSq[j];

                for (int r = 0; r < m; r++) {
                    b_tilde[i][r] -= mu[i][j] * b_tilde[j][r];
                }
            }

            // Calc norm
            double norm = 0;
            for (int r = 0; r < m; r++) norm += b_tilde[i][r] * b_tilde[i][r];
            normsSq[i] = norm;
            mu[i][i] = 1.0;
        }
    }
}
