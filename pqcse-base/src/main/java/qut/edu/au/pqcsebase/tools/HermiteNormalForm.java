package qut.edu.au.pqcsebase.tools;

import java.math.BigInteger;
import java.util.Arrays;

public class HermiteNormalForm {

    /**
     * Compute Hermite Normal Form (HNF).
     * Input and Output are BigInteger[][] to avoid overflow.
     * Core process: Result = Transpose( LowerHNF( Transpose(A) ) )
     */
    public static BigInteger[][] hermiteNormalForm(BigInteger[][] A) {
        if (A.length == 0) return new BigInteger[0][0];

        // 1. Transpose input matrix (A -> A^T)
        BigInteger[][] matrix = transposeBig(A);

        // 2. Compute Lower Triangular HNF
        hnfLower(matrix);

        // 3. Transpose again to get final Upper HNF result
        return transposeBig(matrix);
    }

    /**
     * Compute Lower Triangular HNF.
     * Algorithm: Process columns from right to left (n-1 -> 0).
     */
    private static void hnfLower(BigInteger[][] M) {
        int m = M.length;
        int n = M[0].length;

        for (int c = n - 1; c >= 0; c--) {
            int pivotRow = -1;

            // Find pivot: check current diagonal (c,c) first, then search upwards
            if (c < m && !M[c][c].equals(BigInteger.ZERO)) {
                pivotRow = c;
            } else {
                for (int i = Math.min(c, m - 1); i >= 0; i--) {
                    if (!M[i][c].equals(BigInteger.ZERO)) {
                        pivotRow = i;
                        break;
                    }
                }
            }

            if (pivotRow == -1) continue; // Column is all zeros

            // Swap rows to place pivot at (c, c) if possible
            if (pivotRow != c && c < m) {
                swapRowsBig(M, c, pivotRow);
            }
            int currPivotRow = (c < m) ? c : pivotRow;

            // Eliminate elements above the pivot (Rows 0 .. c-1)
            for (int i = c - 1; i >= 0; i--) {
                if (M[i][c].equals(BigInteger.ZERO)) continue;

                BigInteger u = M[currPivotRow][c];
                BigInteger v = M[i][c];

                BigInteger[] eg = extendedGcdBig(u, v);
                BigInteger g = eg[0];
                BigInteger x = eg[1];
                BigInteger y = eg[2];

                BigInteger uDivG = u.divide(g);
                BigInteger vDivG = v.divide(g);

                BigInteger[] rowP = Arrays.copyOf(M[currPivotRow], n);
                BigInteger[] rowI = M[i];

                for (int k = 0; k < n; k++) {
                    BigInteger valP = rowP[k];
                    BigInteger valI = rowI[k];
                    // New Pivot Row = x*RowP + y*RowI
                    M[currPivotRow][k] = valP.multiply(x).add(valI.multiply(y));
                    // New Row I = -v/g * RowP + u/g * RowI
                    M[i][k] = valP.multiply(vDivG.negate()).add(valI.multiply(uDivG));
                }
            }

            // Ensure pivot is positive
            if (M[currPivotRow][c].signum() < 0) {
                multiplyRowBig(M, currPivotRow, BigInteger.valueOf(-1));
            }

            // Reduce elements below the pivot (Rows c+1 .. m-1)
            for (int i = c + 1; i < m; i++) {
                if (M[i][c].equals(BigInteger.ZERO)) continue;
                BigInteger pivot = M[currPivotRow][c];
                BigInteger val = M[i][c];

                // q = floor(val / pivot)
                BigInteger q;
                if (val.compareTo(BigInteger.ZERO) >= 0) {
                    q = val.divide(pivot);
                } else {
                    BigInteger[] dr = val.divideAndRemainder(pivot);
                    q = dr[0];
                    if (dr[1].signum() != 0) q = q.subtract(BigInteger.ONE);
                }

                if (!q.equals(BigInteger.ZERO)) {
                    rowCombineBig(M, i, BigInteger.ONE, currPivotRow, q.negate());
                }
            }
        }
    }

    // --- Utility Functions ---

    private static BigInteger[][] transposeBig(BigInteger[][] A) {
        int m = A.length;
        int n = A[0].length;
        BigInteger[][] T = new BigInteger[n][m];
        for (int i = 0; i < n; i++) Arrays.fill(T[i], BigInteger.ZERO);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                T[j][i] = A[i][j];
            }
        }
        return T;
    }

    private static void swapRowsBig(BigInteger[][] M, int i, int j) {
        BigInteger[] tmp = M[i];
        M[i] = M[j];
        M[j] = tmp;
    }

    private static void multiplyRowBig(BigInteger[][] M, int row, BigInteger k) {
        for (int j = 0; j < M[row].length; j++) {
            M[row][j] = M[row][j].multiply(k);
        }
    }

    private static void rowCombineBig(BigInteger[][] M, int targetRow, BigInteger aCoef, int sourceRow, BigInteger bCoef) {
        int cols = M[0].length;
        for (int j = 0; j < cols; j++) {
            M[targetRow][j] = aCoef.multiply(M[targetRow][j]).add(bCoef.multiply(M[sourceRow][j]));
        }
    }

    private static BigInteger[] extendedGcdBig(BigInteger a, BigInteger b) {
        BigInteger x = BigInteger.ZERO, y = BigInteger.ONE;
        BigInteger lastX = BigInteger.ONE, lastY = BigInteger.ZERO;
        BigInteger temp;
        while (!b.equals(BigInteger.ZERO)) {
            BigInteger[] dr = a.divideAndRemainder(b);
            BigInteger q = dr[0];
            BigInteger r = dr[1];
            a = b; b = r;
            temp = x; x = lastX.subtract(q.multiply(x)); lastX = temp;
            temp = y; y = lastY.subtract(q.multiply(y)); lastY = temp;
        }
        if (a.signum() < 0) {
            a = a.negate(); lastX = lastX.negate(); lastY = lastY.negate();
        }
        return new BigInteger[]{a, lastX, lastY};
    }
}