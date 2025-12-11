package qut.edu.au.pqcsebase.algorithms;

import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;

import static qut.edu.au.pqcsebase.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebase.tools.BigIntMatrix.multiplyMatrixVector;

public class SampleRight {

    /**
     *
     * @param A a random matrix no trapdoor basis n*k mod q
     * @param B a matrix have trapdoor basis TB   n*m mod q
     * @param R a matrix R k*m mod q  a random matrix in {1, -1} m*m
     * @param TB lattice B trapdoor basis m*m mod q
     * @param u target vector length n  n*1 mod q
     * @param s Gaussian parameter
     * @param nParam security parameter
     * @return e a short vector length m + k such that F2 = [A|AR + B]   F*e = u mod q  (m + k)*1     */
    public static BigInteger[] sampleRight(BigIntMatrix A, BigIntMatrix B, BigIntMatrix R, BigIntMatrix TB, BigInteger[] u, double s, long nParam) {
        int n = A.getRowDimension();
        int k = A.getColumnDimension();
        int m = B.getColumnDimension();
        BigInteger q = A.getModulus();

        //1 construct F2 = [A | AR + B] n*(m+k)
        // compute AR n*m
        BigIntMatrix AR = A.multiplyQ(R);
        // compute AR + B n*m
        BigIntMatrix AR_plus_B = AR.addQ(B);
        // construct F2 n*(m+k)
        BigIntMatrix F2 = BigIntMatrix.columnConcat(A, AR_plus_B);

        //2.Construct the extended basis T_F2 (m+k,m+k)  F2 * T_F2=0 mod q
        BigIntMatrix T_F2 = new BigIntMatrix(m + k, m + k, q);

        //2.1 Map basis vectors from T_B (i = 1 to m)  t_i = [-R * b_i ; b_i]  (m+k)*1
        for (int i = 0; i < m; i++) {
            // get b_i from T_B
            BigInteger[] b_i = new BigInteger[m];
            for (int col = 0; col < m; col++) {
                b_i[col] = TB.get(col, i);
            }
            // compute -R * b_i  k*1
            BigInteger[] topPart = multiplyMatrixVector(R, b_i);
            // Negate it: - (R * b_i)
            for (int j = 0; j < k; j++) {
                topPart[j] = topPart[j].negate().mod(q);
            }

            // Left part (size k): -R * b_i
            for (int j = 0; j < k; j++) T_F2.set(i, j, topPart[j]);
            // Right part (size m): b_i
            for (int j = 0; j < m; j++) T_F2.set(i, k + j, b_i[j]);
        }

        //2.2 Extension vectors (i = 1 to k) t_{m+i} = [w_i - R* u_i; u_i]  (m+k)*1
        // w_i (i-th column of Identity_k)
        for (int i = 0; i < k; i++) {

            // B*u_i = -A*w_i  use samplePre to solve for u_i
            BigInteger[] y = new BigInteger[n];
            for (int row = 0; row < n; row++) {
                y[row] = A.get(row, i).negate().mod(q);
            }
            BigInteger[] u_i = samplePre(B, TB, s, y);

            // Top part: w_i - R * u_i
            BigInteger[] R_u_i = multiplyMatrixVector(R, u_i);
            BigInteger[] topPart = new BigInteger[k];
            for (int j = 0; j < k; j++) {
                BigInteger w_val = (j == i) ? BigInteger.ONE : BigInteger.ZERO;
                // w_i - (R * u_i)
                topPart[j] = w_val.subtract(R_u_i[j]).mod(q);
            }

            // Fill T_F2 row (m + i)
            // Left part (size k): w_i - R * u_i
            for (int j = 0; j < k; j++) T_F2.set(m + i, j, topPart[j]);
            // Right part (size m): u_i
            for (int j = 0; j < m; j++) T_F2.set(m + i, k + j, u_i[j]);
        }

        //3 solution e F2*e = u mod q
        return samplePre(F2, T_F2.transpose(), s, u);
    }

}
