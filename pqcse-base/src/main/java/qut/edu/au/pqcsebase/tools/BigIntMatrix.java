package qut.edu.au.pqcsebase.tools;

import lombok.Getter;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

/**
 * BigInteger Matrix Wrapper supporting modular arithmetic (mod q).
 * Replaces LongMatrix to solve overflow issues with large n and q.
 */
public class BigIntMatrix {
    private final int rows;
    private final int cols;
    private final BigInteger q;
    @Getter
    private final BigInteger[][] data;

    public BigIntMatrix(int rows, int cols, BigInteger q) {
        if (rows <= 0 || cols <= 0) throw new IllegalArgumentException("rows/cols must be > 0");
        this.rows = rows;
        this.cols = cols;
        this.q = q;
        this.data = new BigInteger[rows][cols];
        // Initialize to 0
        for (int i = 0; i < rows; i++) {
            Arrays.fill(this.data[i], BigInteger.ZERO);
        }
    }

    // --- Static Factory Methods ---

    public static BigIntMatrix random(int rows, int cols, BigInteger q) {
        BigIntMatrix m = new BigIntMatrix(rows, cols, q);
        SecureRandom random = new SecureRandom();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                // Generate random number in [0, q-1]
                BigInteger val;
                do {
                    val = new BigInteger(q.bitLength(), random);
                } while (val.compareTo(q) >= 0);
                m.data[i][j] = val;
            }
        }
        return m;
    }

    public static BigIntMatrix zeroMatrix(int rows, int cols, BigInteger q) {
        return new BigIntMatrix(rows, cols, q);
    }

    public static BigIntMatrix identity(int k, BigInteger q) {
        BigIntMatrix m = new BigIntMatrix(k, k, q);
        for (int i = 0; i < k; i++) {
            m.data[i][i] = BigInteger.ONE;
        }
        return m;
    }

    public static BigIntMatrix diagonal(int dimension, BigInteger content, BigInteger q) {
        BigIntMatrix m = new BigIntMatrix(dimension, dimension, q);
        for (int i = 0; i < dimension; i++) {
            m.data[i][i] = content;
        }
        return m;
    }

    public static BigIntMatrix fromArray(long[][] demo, BigInteger bigInteger) {

        int rows = demo.length;
        int cols = demo[0].length;
        BigIntMatrix m = new BigIntMatrix(rows, cols, bigInteger);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                m.data[i][j] = BigInteger.valueOf(demo[i][j]);
            }
        }
        return m;
    }

    // --- Getter / Setter ---

    public int getRowDimension() { return rows; }
    public int getColumnDimension() { return cols; }
    public BigInteger getModulus() { return q; }

    public BigInteger get(int r, int c) {
        return data[r][c];
    }

    public void set(int r, int c, BigInteger value) {
        this.data[r][c] = value;
    }

    public void set(int r, int c, long value) {
        this.data[r][c] = BigInteger.valueOf(value);
    }

    // --- Operations ---

    // Modular addition
    public BigIntMatrix addQ(BigIntMatrix other) {
        checkSameShape(other);
        BigIntMatrix res = new BigIntMatrix(rows, cols, q);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res.data[i][j] = this.data[i][j].add(other.data[i][j]).mod(q);
            }
        }
        return res;
    }

    // Standard addition (no mod q)
    public BigIntMatrix add(BigIntMatrix other) {
        checkSameShape(other);
        BigIntMatrix res = new BigIntMatrix(rows, cols, q);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res.data[i][j] = this.data[i][j].add(other.data[i][j]);
            }
        }
        return res;
    }

    // Modular subtraction
    public BigIntMatrix subtractQ(BigIntMatrix other) {
        checkSameShape(other);
        BigIntMatrix res = new BigIntMatrix(rows, cols, q);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res.data[i][j] = this.data[i][j].subtract(other.data[i][j]).mod(q);
            }
        }
        return res;
    }

    // Standard subtraction (no mod q)
    public BigIntMatrix subtract(BigIntMatrix other) {
        checkSameShape(other);
        BigIntMatrix res = new BigIntMatrix(rows, cols, q);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res.data[i][j] = this.data[i][j].subtract(other.data[i][j]);
            }
        }
        return res;
    }

    // Modular multiplication
    public BigIntMatrix multiplyQ(BigIntMatrix other) {
        if (this.cols != other.rows) throw new IllegalArgumentException("Dimension mismatch");
        if (!this.q.equals(other.q)) throw new IllegalArgumentException("Moduli mismatch");

        BigIntMatrix res = new BigIntMatrix(this.rows, other.cols, q);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                BigInteger sum = BigInteger.ZERO;
                for (int k = 0; k < this.cols; k++) {
                    sum = sum.add(this.data[i][k].multiply(other.data[k][j]));
                }
                res.data[i][j] = sum.mod(q);
            }
        }
        return res;
    }

    // Standard multiplication (no mod q)
    public BigIntMatrix multiply(BigIntMatrix other) {
        if (this.cols != other.rows) throw new IllegalArgumentException("Dimension mismatch");

        BigIntMatrix res = new BigIntMatrix(this.rows, other.cols, q);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                BigInteger sum = BigInteger.ZERO;
                for (int k = 0; k < this.cols; k++) {
                    sum = sum.add(this.data[i][k].multiply(other.data[k][j]));
                }
                res.data[i][j] = sum;
            }
        }
        return res;
    }

    // Transpose
    public BigIntMatrix transpose() {
        BigIntMatrix t = new BigIntMatrix(cols, rows, q);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                t.data[j][i] = this.data[i][j];
            }
        }
        return t;
    }

    // Get submatrix
    public BigIntMatrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) {
        int r = endRow - startRow + 1;
        int c = endCol - startCol + 1;
        BigIntMatrix sub = new BigIntMatrix(r, c, this.q);
        for (int i = 0; i < r; i++) {
            for(int j = 0; j < c; j++) {
                sub.data[i][j] = this.data[startRow + i][startCol + j];
            }
        }
        return sub;
    }

    // --- Static Tools ---

    public static BigIntMatrix rowConcat(BigIntMatrix A, BigIntMatrix B) {
        if (A.cols != B.cols) throw new IllegalArgumentException("Cols must match");
        BigIntMatrix res = new BigIntMatrix(A.rows + B.rows, A.cols, A.q);
        for (int i = 0; i < A.rows; i++) System.arraycopy(A.data[i], 0, res.data[i], 0, A.cols);
        for (int i = 0; i < B.rows; i++) System.arraycopy(B.data[i], 0, res.data[A.rows + i], 0, A.cols);
        return res;
    }

    public static BigIntMatrix columnConcat(BigIntMatrix A, BigIntMatrix B) {
        if (A.rows != B.rows) throw new IllegalArgumentException("Rows must match");
        BigIntMatrix res = new BigIntMatrix(A.rows, A.cols + B.cols, A.q);
        for (int i = 0; i < A.rows; i++) {
            System.arraycopy(A.data[i], 0, res.data[i], 0, A.cols);
            System.arraycopy(B.data[i], 0, res.data[i], A.cols, B.cols);
        }
        return res;
    }

    /**
     * Compute the basis of the image of matrix A (mod q).
     * Uses Gaussian elimination over BigInteger.
     */
    public static BigIntMatrix image(BigIntMatrix A) {
        BigInteger q = A.q;
        int rows = A.rows;
        int cols = A.cols;

        // Deep copy and reduce mod q
        BigInteger[][] mat = new BigInteger[rows][cols];
        for(int i=0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                mat[i][j] = A.data[i][j].mod(q);
            }
        }

        int r = 0;
        for (int c = 0; c < cols && r < rows; c++) {
            int piv = r;
            while (piv < rows && mat[piv][c].equals(BigInteger.ZERO)) piv++;
            if (piv == rows) continue;

            // Swap rows
            if (piv != r) {
                BigInteger[] tmp = mat[r];
                mat[r] = mat[piv];
                mat[piv] = tmp;
            }

            BigInteger pivot = mat[r][c];
            // Compute inverse: pivot^-1 mod q
            BigInteger invPivot;
            try {
                invPivot = pivot.modInverse(q);
            } catch (ArithmeticException e) {
                continue;
            }

            for (int i = r + 1; i < rows; i++) {
                if (mat[i][c].equals(BigInteger.ZERO)) continue;
                BigInteger factor = mat[i][c].multiply(invPivot).mod(q);
                for (int j = c; j < cols; j++) {
                    mat[i][j] = mat[i][j].subtract(factor.multiply(mat[r][j])).mod(q);
                }
            }
            r++;
        }

        // Extract first r rows
        BigIntMatrix X = new BigIntMatrix(r, cols, q);
        for(int i=0; i<r; i++) {
            System.arraycopy(mat[i], 0, X.data[i], 0, cols);
        }
        return X;
    }

    /**
     * Matrix inversion modulo q (Gauss-Jordan elimination).
     * Corresponds to LongMatrix.inverseQ.
     */
    public static BigIntMatrix inverseQ(BigIntMatrix A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("matrix must be square");
        BigInteger q = A.q;

        // Construct augmented matrix n x 2n: [A | I]
        BigInteger[][] aug = new BigInteger[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                aug[i][j] = A.data[i][j].mod(q);
            }
            for (int j = 0; j < n; j++) {
                aug[i][n + j] = (i == j) ? BigInteger.ONE : BigInteger.ZERO;
            }
        }

        for (int i = 0; i < n; i++) {
            // Find an invertible pivot row (starting from i) where gcd(pivot, q) == 1
            int piv = i;
            while (piv < n) {
                BigInteger val = aug[piv][i];
                // Check if coprime with q
                if (!val.gcd(q).equals(BigInteger.ONE)) {
                    piv++;
                } else {
                    break;
                }
            }

            if (piv == n) {
                throw new ArithmeticException("matrix not invertible mod " + q);
            }

            // Swap rows
            if (piv != i) {
                BigInteger[] tmp = aug[i];
                aug[i] = aug[piv];
                aug[piv] = tmp;
            }

            // Normalize pivot row: row[i] = row[i] * pivot^-1 mod q
            BigInteger pivot = aug[i][i];
            BigInteger invPivot = pivot.modInverse(q);
            for (int j = 0; j < 2 * n; j++) {
                aug[i][j] = aug[i][j].multiply(invPivot).mod(q);
            }

            // Eliminate other rows
            for (int r = 0; r < n; r++) {
                if (r == i) continue;
                BigInteger factor = aug[r][i];
                if (factor.equals(BigInteger.ZERO)) continue;

                for (int j = 0; j < 2 * n; j++) {
                    // aug[r][j] = aug[r][j] - factor * aug[i][j]
                    BigInteger val = aug[r][j].subtract(factor.multiply(aug[i][j])).mod(q);
                    aug[r][j] = val;
                }
            }
        }

        // Extract the right half as the inverse matrix
        BigIntMatrix inv = new BigIntMatrix(n, n, q);
        for (int i = 0; i < n; i++) {
            System.arraycopy(aug[i], n, inv.data[i], 0, n);
        }
        return inv;
    }

    /**
     * Exact matrix inversion over integers (no modulo).
     * Only applicable for Unimodular Matrices (det = +/- 1) or specific invertible integer matrices.
     * Corresponds to LongMatrix.inverse.
     */
    public static BigIntMatrix inverse(BigIntMatrix A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("matrix must be square");

        // Construct augmented matrix n x 2n
        BigInteger[][] aug = new BigInteger[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                aug[i][j] = A.data[i][j];
            }
            for (int j = 0; j < n; j++) {
                aug[i][n + j] = (i == j) ? BigInteger.ONE : BigInteger.ZERO;
            }
        }

        for (int i = 0; i < n; i++) {
            // Find non-zero pivot
            int piv = i;
            while (piv < n && aug[piv][i].equals(BigInteger.ZERO)) piv++;

            if (piv == n) throw new ArithmeticException("matrix not invertible (singular)");

            // Swap rows
            if (piv != i) {
                BigInteger[] tmp = aug[i];
                aug[i] = aug[piv];
                aug[piv] = tmp;
            }

            // Normalize pivot row (must be exactly divisible)
            BigInteger pivot = aug[i][i];

            for (int j = 0; j < 2 * n; j++) {
                BigInteger[] divAndRem = aug[i][j].divideAndRemainder(pivot);
                if (!divAndRem[1].equals(BigInteger.ZERO)) {
                    throw new ArithmeticException("matrix not invertible over integers (division not exact)");
                }
                aug[i][j] = divAndRem[0];
            }

            // Eliminate other rows
            for (int r = 0; r < n; r++) {
                if (r == i) continue;
                BigInteger factor = aug[r][i];
                if (factor.equals(BigInteger.ZERO)) continue;

                for (int j = 0; j < 2 * n; j++) {
                    aug[r][j] = aug[r][j].subtract(factor.multiply(aug[i][j]));
                }
            }
        }

        // Extract the right half
        BigIntMatrix inv = new BigIntMatrix(n, n, A.q);
        for (int i = 0; i < n; i++) {
            System.arraycopy(aug[i], n, inv.data[i], 0, n);
        }
        return inv;
    }

    /**
     * Gaussian elimination over integers (no modulo), returns the rank.
     * Uses the strategy: row_i := pivot * row_i - factor * row_r to avoid division.
     * WARNING: Values grow exponentially. BigInteger is essential here.
     * This method modifies the matrix in-place.
     */
    public static int gauss(BigIntMatrix A) {
        int rows = A.rows;
        int cols = A.cols;
        BigInteger[][] data = A.data;

        int r = 0;
        for (int c = 0; c < cols && r < rows; c++) {
            // Find pivot row
            int piv = r;
            while (piv < rows && data[piv][c].equals(BigInteger.ZERO)) piv++;
            if (piv == rows) continue;

            // Swap rows if necessary
            if (piv != r) {
                BigInteger[] tmp = data[r];
                data[r] = data[piv];
                data[piv] = tmp;
            }

            BigInteger pivot = data[r][c];

            // Eliminate rows below
            for (int i = r + 1; i < rows; i++) {
                if (data[i][c].equals(BigInteger.ZERO)) continue;
                BigInteger factor = data[i][c];

                // Update row i: A[i][j] = pivot * A[i][j] - factor * A[r][j]
                for (int j = c; j < cols; j++) {
                    BigInteger term1 = pivot.multiply(data[i][j]);
                    BigInteger term2 = factor.multiply(data[r][j]);
                    data[i][j] = term1.subtract(term2);
                }
            }
            r++;
        }
        return r;
    }


    private void checkSameShape(BigIntMatrix other) {
        if (this.rows != other.rows || this.cols != other.cols) throw new IllegalArgumentException("Shape mismatch");
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (BigInteger[] row : data) {
            sb.append(Arrays.toString(row)).append(System.lineSeparator());
        }
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if(this == o) return true;
        if(o == null || getClass() != o.getClass()) return false;
        BigIntMatrix that = (BigIntMatrix) o;
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                if(!this.data[i][j].mod(q).equals(that.data[i][j].mod(q))) return false;
            }
        }
        return true;
    }

    public String isZeroMatrix() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (!data[i][j].mod(q).equals(BigInteger.ZERO)) {
                    return "Not a zero matrix.";
                }
            }
        }
        return "Is a zero matrix.";
    }
}

