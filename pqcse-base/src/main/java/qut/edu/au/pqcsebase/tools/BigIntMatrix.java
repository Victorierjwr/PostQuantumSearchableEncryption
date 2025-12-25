package qut.edu.au.pqcsebase.tools;

import lombok.Getter;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.stream.IntStream;

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

        // 使用 parallelStream 并行计算每一行
        IntStream.range(0, this.rows).parallel().forEach(i -> {
            for (int j = 0; j < other.cols; j++) {
                BigInteger sum = BigInteger.ZERO;
                for (int k = 0; k < this.cols; k++) {

                    sum = sum.add(this.data[i][k].multiply(other.data[k][j]));
                }
                res.data[i][j] = sum.mod(q);
            }
        });
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
    /** Validate same shape
     */
    public static BigIntMatrix rowConcat(BigIntMatrix A, BigIntMatrix B) {
        if (A.cols != B.cols) throw new IllegalArgumentException("Cols must match");
        BigIntMatrix res = new BigIntMatrix(A.rows + B.rows, A.cols, A.q);
        for (int i = 0; i < A.rows; i++) System.arraycopy(A.data[i], 0, res.data[i], 0, A.cols);
        for (int i = 0; i < B.rows; i++) System.arraycopy(B.data[i], 0, res.data[A.rows + i], 0, A.cols);
        return res;
    }

    /** Column concatenation: [A ; B]
     */
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

    /**
     * Compute the rank of the matrix modulo a prime p.
     * Much faster and avoids OutOfMemoryError compared to integer Gaussian elimination.
     * default BigInteger primeP = BigInteger.valueOf(1000000007);
     */
    public static int gaussMod(BigIntMatrix A,BigInteger primeP) {
        BigInteger p = BigInteger.valueOf(1000000007);
        if (!primeP.equals(BigInteger.ZERO)) p = primeP;
        int rows = A.getRowDimension();
        int cols = A.getColumnDimension();

        BigInteger[][] mat = new BigInteger[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                mat[i][j] = A.get(i, j).mod(p);
            }
        }

        int rank = 0;
        int[] pivotRow = new int[cols];
        java.util.Arrays.fill(pivotRow, -1);

        for (int j = 0; j < cols && rank < rows; j++) {
            int sel = -1;
            for (int i = rank; i < rows; i++) {
                if (!mat[i][j].equals(BigInteger.ZERO)) {
                    sel = i;
                    break;
                }
            }

            if (sel == -1) {
                continue; //all zero column, skip
            }

            // 2. 交换行 (Swap)
            if (sel != rank) {
                BigInteger[] temp = mat[rank];
                mat[rank] = mat[sel];
                mat[sel] = temp;
            }

            BigInteger inv = mat[rank][j].modInverse(p);

            for (int i = 0; i < rows; i++) {
                if (i != rank && !mat[i][j].equals(BigInteger.ZERO)) {
                    // factor = mat[i][j] * inv mod p
                    BigInteger factor = mat[i][j].multiply(inv).mod(p);
                    for (int k = j; k < cols; k++) {
                        // mat[i][k] = mat[i][k] - factor * mat[rank][k]
                        BigInteger sub = factor.multiply(mat[rank][k]).mod(p);
                        mat[i][k] = mat[i][k].subtract(sub).mod(p);
                        if (mat[i][k].signum() < 0) mat[i][k] = mat[i][k].add(p);
                    }
                }
            }
            rank++;
        }
        return rank;
    }

    /**
     * Compute the determinant of matrix A using Bareiss algorithm (no modulo).
     * More efficient and numerically stable than standard LU decomposition for integer matrices.
     */
    public static BigInteger determinant(BigIntMatrix A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("matrix must be square");
        if (n == 0) return BigInteger.ZERO;
        // 深拷贝原矩阵（不取模）
        BigInteger[][] mat = new BigInteger[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A.data[i], 0, mat[i], 0, n);
        }

        int sign = 1;
        if (n == 1) {
            return mat[0][0];
        }

        for (int k = 0; k < n - 1; k++) {
            // 确保枢轴非零（必要时交换行）
            if (mat[k][k].equals(BigInteger.ZERO)) {
                int r = k + 1;
                while (r < n && mat[r][k].equals(BigInteger.ZERO)) r++;
                if (r == n) return BigInteger.ZERO; // 整体奇异
                BigInteger[] tmp = mat[k];
                mat[k] = mat[r];
                mat[r] = tmp;
                sign = -sign;
            }
            BigInteger pivot = mat[k][k];
            BigInteger denom = (k == 0) ? BigInteger.ONE : mat[k - 1][k - 1];

            // Bareiss 更新：保证整除
            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    BigInteger num = mat[i][j].multiply(pivot).subtract(mat[i][k].multiply(mat[k][j]));
                    // 除法在 Bareiss 中应当是整除
                    mat[i][j] = (denom.equals(BigInteger.ONE)) ? num : num.divide(denom);
                }
                mat[i][k] = BigInteger.ZERO;
            }
        }

        BigInteger det = mat[n - 1][n - 1];
        return (sign == -1) ? det.negate() : det;
    }

    /**
     * Compute the determinant of matrix A modulo q using modified Gaussian elimination.
     */
    public static BigInteger determinantMod(BigIntMatrix A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("matrix must be square");
        BigInteger q = A.q;
        if (n == 0) return BigInteger.ZERO;

        // 深拷贝并取模
        BigInteger[][] m = new BigInteger[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m[i][j] = A.data[i][j].mod(q);
            }
        }

        BigInteger det = BigInteger.ONE;
        int sign = 1;

        for (int i = 0; i < n; i++) {
            // 找到与 q 互质的枢轴行
            int piv = i;
            while (piv < n) {
                BigInteger val = m[piv][i];
                if (!val.equals(BigInteger.ZERO) && val.gcd(q).equals(BigInteger.ONE)) break;
                piv++;
            }
            if (piv == n) {
                // 若没有与 q 互质的枢轴，则行列式在模 q 下为 0
                return BigInteger.ZERO;
            }
            if (piv != i) {
                BigInteger[] tmp = m[i];
                m[i] = m[piv];
                m[piv] = tmp;
                sign = -sign;
            }

            BigInteger pivot = m[i][i].mod(q);
            det = det.multiply(pivot).mod(q);

            // 规范化当前行（变为 1）并消去下面行
            BigInteger invPivot = pivot.modInverse(q);
            for (int j = i; j < n; j++) {
                m[i][j] = m[i][j].multiply(invPivot).mod(q);
            }
            for (int r = i + 1; r < n; r++) {
                BigInteger factor = m[r][i];
                if (factor.equals(BigInteger.ZERO)) continue;
                for (int c = i; c < n; c++) {
                    m[r][c] = m[r][c].subtract(factor.multiply(m[i][c])).mod(q);
                    if (m[r][c].signum() < 0) m[r][c] = m[r][c].add(q);
                }
            }
        }

        if (sign == -1) det = det.negate().mod(q);
        return det.mod(q);
    }

    /**
     * verfy ||S|| <= L
     * Compute the GSN length of the basis represented by the columns of the matrix.
     * Uses the Euclidean norm of each column and returns the maximum.
     */
    public static BigInteger GSN_Length(BigIntMatrix matrix) {
        BigInteger maxLen = BigInteger.ZERO;

        for (int j = 0; j < matrix.getColumnDimension(); j++) {
            // 使用 BigInteger sum of squares
            BigInteger sumSq = BigInteger.ZERO;
            for (int i = 0; i < matrix.getRowDimension(); i++) {
                BigInteger val = matrix.get(i, j);
                sumSq = sumSq.add(val.multiply(val));
            }

            BigInteger colNorm = sqrt(sumSq);

            if (colNorm.compareTo(maxLen) > 0) {
                maxLen = colNorm;
            }
        }
        return maxLen;
    }

    /**
     * Compute the integer square root of a BigInteger n.
     * Uses binary search method.
     */
    private static BigInteger sqrt(BigInteger n) {
        if (n.signum() < 0) throw new IllegalArgumentException("Negative argument.");
        if (n.equals(BigInteger.ZERO) || n.equals(BigInteger.ONE)) return n;

        BigInteger a = BigInteger.ONE;
        BigInteger b = n.shiftRight(5).add(BigInteger.valueOf(8));
        while (b.compareTo(a) >= 0) {
            BigInteger mid = a.add(b).shiftRight(1);
            if (mid.multiply(mid).compareTo(n) > 0) b = mid.subtract(BigInteger.ONE);
            else a = mid.add(BigInteger.ONE);
        }
        return a.subtract(BigInteger.ONE);
    }

    /**
     * Helper: Multiply Matrix (rows x cols) by Vector (cols).
     * Returns a vector of length 'rows'.
     */
    public static BigInteger[] multiplyMatrixVector(BigIntMatrix M, BigInteger[] v) {
        int rows = M.getRowDimension();
        int cols = M.getColumnDimension();
        BigInteger q = M.getModulus();

        if (v.length != cols) throw new IllegalArgumentException("Dimension mismatch");

        BigInteger[] result = new BigInteger[rows];
        for (int i = 0; i < rows; i++) {
            BigInteger sum = BigInteger.ZERO;
            for (int j = 0; j < cols; j++) {
                sum = sum.add(M.get(i, j).multiply(v[j]));
            }
            result[i] = sum.mod(q);
        }
        return result;
    }

    /**
     * Swap two columns at indices i and j.
     * @param i Index of the first column.
     * @param j Index of the second column.
     */
    public void swapCols(int i, int j) {
        for (int k = 0; k < this.rows; k++) {
            BigInteger temp = this.data[k][i];
            this.data[k][i] = this.data[k][j];
            this.data[k][j] = temp;
        }
    }

    /**
     * Subtract a multiple of column j from column i.
     * Operation: Col_i = Col_i - factor * Col_j
     *
     * @param i      Target column index.
     * @param j      Source column index.
     * @param factor The integer multiplier.
     */
    public void colSub(int i, int j, BigInteger factor) {
        for (int k = 0; k < this.rows; k++) {
            BigInteger val = this.data[k][j].multiply(factor);
            this.data[k][i] = this.data[k][i].subtract(val);
        }
    }

    /**
     * Creates a deep copy of the matrix.
     * @return A new BigIntMatrix instance with copied data.
     */
    public BigIntMatrix copy() {
        BigIntMatrix newMat = new BigIntMatrix(this.rows, this.cols, this.q);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                newMat.set(i, j, this.data[i][j]);
            }
        }
        return newMat;
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

