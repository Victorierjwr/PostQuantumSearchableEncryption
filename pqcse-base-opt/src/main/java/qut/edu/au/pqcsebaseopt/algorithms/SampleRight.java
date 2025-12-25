package qut.edu.au.pqcsebaseopt.algorithms;

import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.stream.IntStream;

import static qut.edu.au.pqcsebaseopt.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebaseopt.tools.BigIntMatrix.multiplyMatrixVector;

public class SampleRight {

    public static BigInteger[] sampleRight(BigIntMatrix A, BigIntMatrix B, BigIntMatrix R, BigIntMatrix TB, BigInteger[] u, double sigma, long nParam) {
        int n = A.getRows();
        int k = A.getCols();
        int m = B.getCols();
        BigInteger q = A.getModulus();

        // 1. Construct F2 = [A | AR + B]
        // Optimized matrix arithmetic handles this efficiently
        BigIntMatrix AR = A.multiplyQ(R);
        BigIntMatrix AR_plus_B = AR.addQ(B);
        BigIntMatrix F2 = BigIntMatrix.columnConcat(A, AR_plus_B);

        // 2. Construct Extended Basis T_F2
        BigIntMatrix T_F2 = new BigIntMatrix(m + k, m + k, q);

        // 2.1 First m rows: t_i = [-R * b_i ; b_i]
        // Parallelize over m columns of TB
        IntStream.range(0, m).parallel().forEach(i -> {
            // Extract column b_i
            BigInteger[] b_i = new BigInteger[m];
            for (int col = 0; col < m; col++) {
                b_i[col] = TB.get(col, i);
            }
            // Compute -R * b_i
            BigInteger[] topPart = multiplyMatrixVector(R, b_i);
            // Negate and fill
            for (int j = 0; j < k; j++) {
                T_F2.set(i, j, topPart[j].negate().mod(q));
            }
            for (int j = 0; j < m; j++) {
                T_F2.set(i, k + j, b_i[j]);
            }
        });

        // 2.2 Extension vectors (i = 0 to k-1)
        // Parallelize the heavy SamplePre calls!
        IntStream.range(0, k).parallel().forEach(i -> {
            // Solve B*u_i = -A*w_i (where w_i is i-th col of I)
            // -A*w_i is simply the negated i-th column of A
            BigInteger[] y = new BigInteger[n];
            for (int row = 0; row < n; row++) {
                y[row] = A.get(row, i).negate().mod(q);
            }

            // Expensive call, running in parallel
            BigInteger[] u_i = samplePre(B, TB, sigma, y);

            // Compute w_i - R * u_i
            BigInteger[] R_u_i = multiplyMatrixVector(R, u_i);

            // Fill row m+i
            // Left part: w_i - R*u_i
            for (int j = 0; j < k; j++) {
                BigInteger w_val = (j == i) ? BigInteger.ONE : BigInteger.ZERO;
                BigInteger val = w_val.subtract(R_u_i[j]).mod(q);
                T_F2.set(m + i, j, val);
            }
            // Right part: u_i
            for (int j = 0; j < m; j++) {
                T_F2.set(m + i, k + j, u_i[j]);
            }
        });

        // 3. Solve F2 * e = u
        return samplePre(F2, T_F2.transpose(), sigma, u);
    }

    public static BigInteger[] sampleRightFast(BigIntMatrix A, BigIntMatrix B, BigIntMatrix R, BigIntMatrix TB, BigInteger[] u, double sigma, long nParam) {
        int n = A.getRows();
        int k = A.getCols();
        int m = B.getCols();
        BigInteger q = A.getModulus();

        // 1. Construct F2 = [A | AR + B]
        BigIntMatrix AR = A.multiplyQ(R);
        BigIntMatrix AR_plus_B = AR.addQ(B);
        BigIntMatrix F2 = BigIntMatrix.columnConcat(A, AR_plus_B);

        // 2. Construct Extended Basis T_F2
        BigIntMatrix T_F2 = new BigIntMatrix(m + k, m + k, q);

        // 2.1 First m rows
        IntStream.range(0, m).parallel().forEach(i -> {
            BigInteger[] b_i = new BigInteger[m];
            for (int col = 0; col < m; col++) b_i[col] = TB.get(col, i);
            BigInteger[] topPart = multiplyMatrixVector(R, b_i);
            for (int j = 0; j < k; j++) T_F2.set(i, j, topPart[j].negate().mod(q));
            for (int j = 0; j < m; j++) T_F2.set(i, k + j, b_i[j]);
        });

        // =========================================================
        // [CRITICAL OPTIMIZATION] PRE-COMPUTE GSO for TB to avoid OOM
        // =========================================================
        int tb_cols = TB.getCols();
        int tb_rows = TB.getRows();
        double[][] TB_double = new double[tb_cols][tb_rows];

        // Convert to double in parallel
        IntStream.range(0, tb_cols).parallel().forEach(i -> {
            for (int j = 0; j < tb_rows; j++) TB_double[i][j] = TB.get(j, i).doubleValue();
        });

        double[][] TB_tilde = new double[tb_cols][tb_rows];
        double[] normsSq = new double[tb_cols];

        // Compute GSO ONCE
        qut.edu.au.pqcsebaseopt.algorithms.SampleD.computeGramSchmidt(TB_double, TB_tilde, normsSq);
        // =========================================================

        // 2.2 Extension vectors (using precomputed GSO)
        IntStream.range(0, k).parallel().forEach(i -> {
            BigInteger[] y = new BigInteger[n];
            for (int row = 0; row < n; row++) {
                y[row] = A.get(row, i).negate().mod(q);
            }

            // USE FAST VERSION
            BigInteger[] u_i = qut.edu.au.pqcsebaseopt.algorithms.SamplePre.samplePreFast(
                    B, TB, TB_tilde, normsSq, TB_double, sigma, y
            );

            BigInteger[] R_u_i = multiplyMatrixVector(R, u_i);

            for (int j = 0; j < k; j++) {
                BigInteger w_val = (j == i) ? BigInteger.ONE : BigInteger.ZERO;
                BigInteger val = w_val.subtract(R_u_i[j]).mod(q);
                T_F2.set(m + i, j, val);
            }
            for (int j = 0; j < m; j++) {
                T_F2.set(m + i, k + j, u_i[j]);
            }
        });

        // 3. Solve F2 * e = u
        // (This runs once, so standard samplePre is fine, or optimize similarly if needed)
        return samplePre(F2, T_F2.transpose(), sigma, u);
    }
}