package qut.edu.au.pqcsebaseopt.tools.back;

import lombok.Getter;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * High-Performance Hybrid BigInteger Matrix.
 * * Optimization Highlights:
 * 1. Storage: Flattened 1D BigInteger[] array for cache locality and huge number support.
 * 2. Multiplication: 'multiplyQ' uses hybrid approach (BigInt -> long[] -> AVX-like intrinsics) for speed.
 * 3. Compatibility: Retains ALL original methods (gauss, determinant, inverse, etc.) adapted for 1D storage.
 */
public class BigIntMatrix1 {
    @Getter
    private final int rows;
    @Getter
    private final int cols;

    // Primary storage: 1D Flattened BigInteger array
    // Index = row * cols + col
    private final BigInteger[] data;

    @Getter
    private final BigInteger modulus; // Original 'q' (renamed from 'q' to 'modulus' for clarity, getter kept as getModulus or q)
    private final long modLong;       // Cached 'q' as primitive long for fast ops

    // --- Constructors ---

    public BigIntMatrix1(int rows, int cols, BigInteger q) {
        if (rows <= 0 || cols <= 0) throw new IllegalArgumentException("rows/cols must be > 0");
        this.rows = rows;
        this.cols = cols;
        this.modulus = q;

        // Cache modulus for fast conversions
        if (q.bitLength() > 63) {
            this.modLong = Long.MAX_VALUE; // Sentinel: too big for long optimization
        } else {
            this.modLong = q.longValueExact();
        }

        this.data = new BigInteger[rows * cols];
        Arrays.fill(this.data, BigInteger.ZERO);
    }

    // --- Static Factory Methods ---

    public static BigIntMatrix1 random(int rows, int cols, BigInteger q) {
        BigIntMatrix1 m = new BigIntMatrix1(rows, cols, q);
        SecureRandom random = new SecureRandom();
        // Parallel random generation
        IntStream.range(0, rows * cols).parallel().forEach(i -> {
            BigInteger val;
            do {
                val = new BigInteger(q.bitLength(), random);
            } while (val.compareTo(q) >= 0);
            m.data[i] = val;
        });
        return m;
    }

    public static BigIntMatrix1 zeroMatrix(int rows, int cols, BigInteger q) {
        return new BigIntMatrix1(rows, cols, q);
    }

    public static BigIntMatrix1 identity(int k, BigInteger q) {
        BigIntMatrix1 m = new BigIntMatrix1(k, k, q);
        for (int i = 0; i < k; i++) {
            m.data[i * k + i] = BigInteger.ONE;
        }
        return m;
    }

    public static BigIntMatrix1 diagonal(int dimension, BigInteger content, BigInteger q) {
        BigIntMatrix1 m = new BigIntMatrix1(dimension, dimension, q);
        for (int i = 0; i < dimension; i++) {
            m.data[i * dimension + i] = content;
        }
        return m;
    }

    public static BigIntMatrix1 fromArray(long[][] demo, BigInteger bigInteger) {
        int rows = demo.length;
        int cols = demo[0].length;
        BigIntMatrix1 m = new BigIntMatrix1(rows, cols, bigInteger);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                m.data[i * cols + j] = BigInteger.valueOf(demo[i][j]);
            }
        }
        return m;
    }

    // --- Getters / Setters ---

    public int getRowDimension() { return rows; }
    public int getColumnDimension() { return cols; }
    public BigInteger getModulus() { return modulus; }

    // Compatibility: return 2D array copy
    public BigInteger[][] getData() {
        BigInteger[][] out = new BigInteger[rows][cols];
        for(int i=0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                out[i][j] = data[i * cols + j];
            }
        }
        return out;
    }

    public BigInteger get(int r, int c) {
        return data[r * cols + c];
    }

    public void set(int r, int c, BigInteger value) {
        this.data[r * cols + c] = value;
    }

    public void set(int r, int c, long value) {
        this.data[r * cols + c] = BigInteger.valueOf(value);
    }

    // Fast primitive getter for optimization
    public long getLong(int r, int c) {
        return data[r * cols + c].longValue();
    }

    // --- Core Operations (Optimized) ---

    // Modular Addition (Parallel)
    public BigIntMatrix1 addQ(BigIntMatrix1 other) {
        checkSameShape(other);
        BigIntMatrix1 res = new BigIntMatrix1(rows, cols, modulus);
        BigInteger m = this.modulus;
        IntStream.range(0, data.length).parallel().forEach(i -> {
            res.data[i] = this.data[i].add(other.data[i]).mod(m);
        });
        return res;
    }

    // Standard Addition (Parallel)
    public BigIntMatrix1 add(BigIntMatrix1 other) {
        checkSameShape(other);
        BigIntMatrix1 res = new BigIntMatrix1(rows, cols, modulus);
        IntStream.range(0, data.length).parallel().forEach(i -> {
            res.data[i] = this.data[i].add(other.data[i]);
        });
        return res;
    }

    // Modular Subtraction (Parallel)
    public BigIntMatrix1 subtractQ(BigIntMatrix1 other) {
        checkSameShape(other);
        BigIntMatrix1 res = new BigIntMatrix1(rows, cols, modulus);
        BigInteger m = this.modulus;
        IntStream.range(0, data.length).parallel().forEach(i -> {
            res.data[i] = this.data[i].subtract(other.data[i]).mod(m);
        });
        return res;
    }

    // Standard Subtraction (Parallel)
    public BigIntMatrix1 subtract(BigIntMatrix1 other) {
        checkSameShape(other);
        BigIntMatrix1 res = new BigIntMatrix1(rows, cols, modulus);
        IntStream.range(0, data.length).parallel().forEach(i -> {
            res.data[i] = this.data[i].subtract(other.data[i]);
        });
        return res;
    }

    /**
     * Optimized Modular Multiplication.
     * Uses hybrid approach: Converts to long[] for speed, then back to BigInteger.
     * FIX: explicit .mod(modBI) before .longValue() to handle inputs larger than 2^63 (like G and S).
     */
    public BigIntMatrix1 multiplyQ(BigIntMatrix1 other) {
        if (this.cols != other.rows) throw new IllegalArgumentException("Dimension mismatch");
        if (!this.modulus.equals(other.modulus)) throw new IllegalArgumentException("Moduli mismatch");

        BigIntMatrix1 res = new BigIntMatrix1(this.rows, other.cols, this.modulus);

        int n = this.cols;
        int m = other.cols;
        int r = this.rows;
        BigInteger modBI = this.modulus;

        // 1. Snapshot Conversion: BigInteger[] -> long[] for speed
        // [CRITICAL FIX]: Must apply .mod(modBI) first because inputs (like S or G) can be HUGE integers.
        // Direct .longValue() would truncate them, leading to wrong results mod q.
        long[] A_flat = new long[this.data.length];
        long[] B_flat = new long[other.data.length];

        IntStream.range(0, this.data.length).parallel().forEach(i ->
                A_flat[i] = this.data[i].mod(modBI).longValue()
        );
        IntStream.range(0, other.data.length).parallel().forEach(i ->
                B_flat[i] = other.data[i].mod(modBI).longValue()
        );

        BigInteger[] C = res.data;

        // 2. High-Speed Multiplication (Parallel + Intrinsics)
        IntStream.range(0, r).parallel().forEach(i -> {
            int rowHeadA = i * n;
            int rowHeadC = i * m;

            for (int j = 0; j < m; j++) {
                long accLo = 0, accMid = 0, accHi = 0;

                for (int k = 0; k < n; k++) {
                    long aVal = A_flat[rowHeadA + k];
                    long bVal = B_flat[k * m + j];

                    long prodHi = Math.multiplyHigh(aVal, bVal);
                    long prodLo = aVal * bVal;

                    if (prodHi < 0) accHi--; // Sign extension correction

                    long prevLo = accLo;
                    accLo += prodLo;
                    boolean carryLo = Long.compareUnsigned(accLo, prevLo) < 0;

                    long prevMid = accMid;
                    accMid += prodHi;
                    if (carryLo) {
                        accMid++;
                        if (accMid == 0) accHi++;
                    }
                    if (Long.compareUnsigned(accMid, prevMid) < 0) {
                        accHi++;
                    }
                }

                // 3. Reduce and store back as BigInteger
                C[rowHeadC + j] = fastReduce192ToBigInt(accHi, accMid, accLo, modBI);
            }
        });
        return res;
    }

    /**
     * Standard Integer Multiplication (Non-Modular).
     * Necessary for TrapGen S-matrix construction logic.
     * Uses BigInteger accumulation to avoid overflow.
     */
    public BigIntMatrix1 multiply(BigIntMatrix1 other) {
        if (this.cols != other.rows) throw new IllegalArgumentException("Dimension mismatch");

        BigIntMatrix1 res = new BigIntMatrix1(this.rows, other.cols, this.modulus);
        int n = this.cols;
        int m = other.cols;
        int r = this.rows;

        BigInteger[] A = this.data;
        BigInteger[] B = other.data;
        BigInteger[] C = res.data;

        IntStream.range(0, r).parallel().forEach(i -> {
            int rowHeadA = i * n;
            int rowHeadC = i * m;

            for (int j = 0; j < m; j++) {
                BigInteger sum = BigInteger.ZERO;
                for (int k = 0; k < n; k++) {
                    sum = sum.add(A[rowHeadA + k].multiply(B[k * m + j]));
                }
                C[rowHeadC + j] = sum;
            }
        });
        return res;
    }

    // Transpose
    public BigIntMatrix1 transpose() {
        BigIntMatrix1 t = new BigIntMatrix1(cols, rows, modulus);
        IntStream.range(0, rows).parallel().forEach(i -> {
            for (int j = 0; j < cols; j++) {
                t.data[j * rows + i] = this.data[i * cols + j];
            }
        });
        return t;
    }

    // Submatrix
    public BigIntMatrix1 getSubMatrix(int startRow, int endRow, int startCol, int endCol) {
        int r = endRow - startRow + 1;
        int c = endCol - startCol + 1;
        BigIntMatrix1 sub = new BigIntMatrix1(r, c, modulus);
        IntStream.range(0, r).parallel().forEach(i -> {
            System.arraycopy(this.data, (startRow + i) * this.cols + startCol, sub.data, i * c, c);
        });
        return sub;
    }

    // --- Complex Algorithms (Adapted for 1D Storage) ---

    /**
     * Compute the basis of the image of matrix A (mod q).
     * Uses Gaussian elimination over BigInteger.
     */
    public static BigIntMatrix1 image(BigIntMatrix1 A) {
        BigInteger q = A.modulus;
        int rows = A.rows;
        int cols = A.cols;

        // Deep copy
        BigIntMatrix1 X = A.copy();
        BigInteger[] mat = X.data;

        int r = 0;
        for (int c = 0; c < cols && r < rows; c++) {
            int piv = r;
            while (piv < rows && mat[piv * cols + c].equals(BigInteger.ZERO)) piv++;
            if (piv == rows) continue;

            if (piv != r) X.swapRows(r, piv);

            BigInteger pivot = mat[r * cols + c];
            BigInteger invPivot;
            try {
                invPivot = pivot.modInverse(q);
            } catch (ArithmeticException e) {
                continue;
            }

            for (int i = r + 1; i < rows; i++) {
                if (mat[i * cols + c].equals(BigInteger.ZERO)) continue;
                BigInteger factor = mat[i * cols + c].multiply(invPivot).mod(q);
                for (int j = c; j < cols; j++) {
                    int idxI = i * cols + j;
                    int idxR = r * cols + j;
                    mat[idxI] = mat[idxI].subtract(factor.multiply(mat[idxR])).mod(q);
                }
            }
            r++;
        }
        return X.getSubMatrix(0, r - 1, 0, cols - 1);
    }

    /**
     * Matrix inversion modulo q (Gauss-Jordan).
     */
    /**
     * Optimized Inverse Mod Q using Parallelism and Primitive Longs.
     * Complexity: O(n^3) but highly optimized for modern CPUs.
     */
    public static BigIntMatrix1 inverseQ(BigIntMatrix1 A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("matrix must be square");

        // 检查模数是否适合 long 优化
        if (A.modulus.bitLength() > 62) {
            // 如果模数太大，只能回退到慢速版本（或者报错）
            // 这里为了演示，假设必须优化，直接抛出建议
            throw new IllegalArgumentException("Modulus too large for long optimization. Use q < 2^62.");
        }

        long mod = A.modulus.longValue();

        // 1. 构造增广矩阵 [A | I] 使用一维 long 数组
        // 尺寸: n * (2n)
        int rowWidth = 2 * n;
        long[] aug = new long[n * rowWidth];

        // 并行初始化增广矩阵
        IntStream.range(0, n).parallel().forEach(i -> {
            int rowOffset = i * rowWidth;
            // 复制 A 的数据 (A 是 n*n)
            // 假设 A.data 是一维数组 [row * n + col]
            // 如果 A.data 是 BigInteger[]，需要转换
            for (int j = 0; j < n; j++) {
                // 注意：这里调用的是 A.get(i, j) 或者 A.data[i*n+j]
                // 建议使用 getLong 避免对象创建
                aug[rowOffset + j] = A.data[i * n + j].longValue();
            }
            // 设置单位矩阵 I 的部分
            // i-th row, (n + i)-th column is 1
            aug[rowOffset + n + i] = 1L;
        });

        // 2. 高斯-若尔当消元 (Gauss-Jordan)
        for (int i = 0; i < n; i++) {

            if (i % 100 == 0) System.out.println("InverseQ optimized loop: " + i + " / " + n);

            int piv = i;
            int rowOffsetPivot = i * rowWidth;

            // 寻找主元 (Pivot Search) - 这一步必须串行
            // 优化：优先找非零元即可 (素数域下非零即在大数域互质，除非 q 不是素数)
            while (piv < n) {
                long val = aug[piv * rowWidth + i];
                // 简单检查非零即可，如果 q 是素数
                // 严格检查 gcd(val, mod) == 1
                if (val != 0 && gcd(val, mod) == 1) break;
                piv++;
            }

            if (piv == n) throw new ArithmeticException("Matrix singular or pivot not coprime to modulus");

            // 交换行 (Swap Rows)
            if (piv != i) {
                int rowOffsetSwap = piv * rowWidth;
                // 并行交换行数据
                IntStream.range(0, rowWidth).parallel().forEach(col -> {
                    long tmp = aug[rowOffsetPivot + col];
                    aug[rowOffsetPivot + col] = aug[rowOffsetSwap + col];
                    aug[rowOffsetSwap + col] = tmp;
                });
            }

            // 归一化主元行 (Normalize Pivot Row)
            long pivotVal = aug[rowOffsetPivot + i];
            long invPivot = modInverse(pivotVal, mod);

            // 并行处理当前行的归一化
            // 优化：只需要从列 i 开始处理，前面的都是 0
            IntStream.range(i, rowWidth).parallel().forEach(j -> {
                aug[rowOffsetPivot + j] = (aug[rowOffsetPivot + j] * invPivot) % mod;
            });

            // 消元其他行 (Eliminate Other Rows)
            // 这是 O(n^3) 的核心部分，必须并行化！
            final int currentI = i;
            IntStream.range(0, n).parallel().forEach(r -> {
                if (r != currentI) {
                    int rowOffsetTarget = r * rowWidth;
                    long factor = aug[rowOffsetTarget + currentI];

                    if (factor != 0) {
                        // row[r] = row[r] - factor * row[i]
                        // 优化：从列 i 开始，因为 i 之前的列主元行都是 0，不会影响结果
                        for (int j = currentI; j < rowWidth; j++) {
                            long sub = (factor * aug[rowOffsetPivot + j]) % mod;
                            long res = aug[rowOffsetTarget + j] - sub;
                            if (res < 0) res += mod;
                            aug[rowOffsetTarget + j] = res;
                        }
                    }
                }
            });
        }

        // 3. 提取结果矩阵 (右半部分)
        BigIntMatrix1 inv = new BigIntMatrix1(n, n, A.modulus);
        // 并行提取
        IntStream.range(0, n).parallel().forEach(i -> {
            int rowOffsetAug = i * rowWidth;
            int rowOffsetRes = i * n;
            for (int j = 0; j < n; j++) {
                // aug 的第 n+j 列是结果的第 j 列
                inv.data[rowOffsetRes + j] = BigInteger.valueOf(aug[rowOffsetAug + n + j]);
            }
        });

        return inv;
    }

    // --- 必要的辅助方法 (Primitive Long) ---

    private static long modInverse(long a, long m) {
        long m0 = m, t, q;
        long x0 = 0, x1 = 1;
        if (m == 1) return 0;
        while (a > 1) {
            if (m == 0) throw new ArithmeticException("Inverse does not exist");
            q = a / m;
            t = m;
            m = a % m;
            a = t;
            t = x0;
            x0 = x1 - q * x0;
            x1 = t;
        }
        if (x1 < 0) x1 += m0;
        return x1;
    }

    private static long gcd(long a, long b) {
        return b == 0 ? a : gcd(b, a % b);
    }

    /**
     * Exact matrix inversion over integers.
     */
    public static BigIntMatrix1 inverse(BigIntMatrix1 A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("matrix must be square");

        // Similar logic to inverseQ but integer division.
        // Simplified adaption for 1D array.
        BigInteger[] aug = new BigInteger[n * 2 * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A.data, i * n, aug, i * 2 * n, n);
            for (int j = 0; j < n; j++) aug[i * 2 * n + n + j] = (i == j) ? BigInteger.ONE : BigInteger.ZERO;
        }

        for (int i = 0; i < n; i++) {
            int piv = i;
            while (piv < n && aug[piv * 2 * n + i].equals(BigInteger.ZERO)) piv++;
            if (piv == n) throw new ArithmeticException("matrix not invertible (singular)");

            if (piv != i) {
                BigInteger[] tmp = new BigInteger[2 * n];
                System.arraycopy(aug, i * 2 * n, tmp, 0, 2 * n);
                System.arraycopy(aug, piv * 2 * n, aug, i * 2 * n, 2 * n);
                System.arraycopy(tmp, 0, aug, piv * 2 * n, 2 * n);
            }

            BigInteger pivot = aug[i * 2 * n + i];
            for (int j = 0; j < 2 * n; j++) {
                BigInteger[] dr = aug[i * 2 * n + j].divideAndRemainder(pivot);
                if (!dr[1].equals(BigInteger.ZERO)) throw new ArithmeticException("Not invertible over integers");
                aug[i * 2 * n + j] = dr[0];
            }

            for (int r = 0; r < n; r++) {
                if (r == i) continue;
                BigInteger factor = aug[r * 2 * n + i];
                if (factor.equals(BigInteger.ZERO)) continue;
                for (int j = 0; j < 2 * n; j++) {
                    aug[r * 2 * n + j] = aug[r * 2 * n + j].subtract(factor.multiply(aug[i * 2 * n + j]));
                }
            }
        }

        BigIntMatrix1 inv = new BigIntMatrix1(n, n, A.modulus);
        for (int i = 0; i < n; i++) System.arraycopy(aug, i * 2 * n + n, inv.data, i * n, n);
        return inv;
    }

    /**
     * Gaussian elimination over integers (returns rank).
     */
    public static int gauss(BigIntMatrix1 A) {
        int rows = A.rows;
        int cols = A.cols;
        BigInteger[] data = A.data;

        int r = 0;
        for (int c = 0; c < cols && r < rows; c++) {
            int piv = r;
            while (piv < rows && data[piv * cols + c].equals(BigInteger.ZERO)) piv++;
            if (piv == rows) continue;

            if (piv != r) A.swapRows(r, piv);

            BigInteger pivot = data[r * cols + c];

            for (int i = r + 1; i < rows; i++) {
                if (data[i * cols + c].equals(BigInteger.ZERO)) continue;
                BigInteger factor = data[i * cols + c];
                for (int j = c; j < cols; j++) {
                    // data[i][j] = pivot * data[i][j] - factor * data[r][j]
                    BigInteger term1 = pivot.multiply(data[i * cols + j]);
                    BigInteger term2 = factor.multiply(data[r * cols + j]);
                    data[i * cols + j] = term1.subtract(term2);
                }
            }
            r++;
        }
        return r;
    }

    /**
     * Rank modulo p.
     */
    public static int gaussMod(BigIntMatrix1 A, BigInteger primeP) {
        BigInteger p;
        if (!primeP.equals(BigInteger.ZERO)) p = primeP;
        else {
            p = BigInteger.valueOf(1000000007);
        }
        int rows = A.rows;
        int cols = A.cols;

        // Local copy mod p
        BigInteger[] mat = new BigInteger[rows * cols];
        IntStream.range(0, rows * cols).parallel().forEach(i -> mat[i] = A.data[i].mod(p));

        int rank = 0;
        for (int j = 0; j < cols && rank < rows; j++) {
            int sel = -1;
            for (int i = rank; i < rows; i++) {
                if (!mat[i * cols + j].equals(BigInteger.ZERO)) {
                    sel = i;
                    break;
                }
            }
            if (sel == -1) continue;

            if (sel != rank) {
                // swap rows
                BigInteger[] temp = new BigInteger[cols];
                System.arraycopy(mat, rank * cols, temp, 0, cols);
                System.arraycopy(mat, sel * cols, mat, rank * cols, cols);
                System.arraycopy(temp, 0, mat, sel * cols, cols);
            }

            BigInteger inv = mat[rank * cols + j].modInverse(p);

            for (int i = 0; i < rows; i++) {
                if (i != rank && !mat[i * cols + j].equals(BigInteger.ZERO)) {
                    BigInteger factor = mat[i * cols + j].multiply(inv).mod(p);
                    for (int k = j; k < cols; k++) {
                        BigInteger sub = factor.multiply(mat[rank * cols + k]).mod(p);
                        mat[i * cols + k] = mat[i * cols + k].subtract(sub).mod(p);
                        if (mat[i * cols + k].signum() < 0) mat[i * cols + k] = mat[i * cols + k].add(p);
                    }
                }
            }
            rank++;
        }
        return rank;
    }

    /**
     * Determinant (Bareiss).
     */
    public static BigInteger determinant(BigIntMatrix1 A) {
        int n = A.rows;
        if (n != A.cols) throw new IllegalArgumentException("Square only");
        if (n == 0) return BigInteger.ZERO;

        // Deep copy
        BigInteger[] mat = Arrays.copyOf(A.data, A.data.length);
        int sign = 1;

        for (int k = 0; k < n - 1; k++) {
            if (mat[k * n + k].equals(BigInteger.ZERO)) {
                int r = k + 1;
                while (r < n && mat[r * n + k].equals(BigInteger.ZERO)) r++;
                if (r == n) return BigInteger.ZERO;
                // Swap
                for(int x=0; x<n; x++) {
                    BigInteger tmp = mat[k*n+x];
                    mat[k*n+x] = mat[r*n+x];
                    mat[r*n+x] = tmp;
                }
                sign = -sign;
            }
            BigInteger pivot = mat[k * n + k];
            BigInteger denom = (k == 0) ? BigInteger.ONE : mat[(k - 1) * n + (k - 1)];

            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    BigInteger num = mat[i*n+j].multiply(pivot).subtract(mat[i*n+k].multiply(mat[k*n+j]));
                    mat[i*n+j] = (denom.equals(BigInteger.ONE)) ? num : num.divide(denom);
                }
                mat[i * n + k] = BigInteger.ZERO;
            }
        }
        BigInteger det = mat[(n - 1) * n + (n - 1)];
        return (sign == -1) ? det.negate() : det;
    }

    public static BigInteger determinantMod(BigIntMatrix1 A) {
        // 如果模数 q 小于 63 bit (Long.MAX_VALUE)，走高速通道
        if (A.getModulus().bitLength() < 63) {
            return determinantModFast(A);
        }
        // 否则走原本的通用通道
        return determinantModGeneral(A);
    }

    // === 高速通道 (基于 long) ===
    private static BigInteger determinantModFast(BigIntMatrix1 A) {
        System.out.println("Using fast determinant mod long channel.");
        int n = A.getRowDimension();
        if (n != A.getColumnDimension()) return BigInteger.ZERO;

        long q = A.getModulus().longValue();

        // 1. 转存为 long 数组 (扁平化一维数组性能更好)
        long[] m = new long[n * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m[i * n + j] = A.get(i, j).longValue();
            }
        }

        long det = 1;
        int sign = 1;

        for (int i = 0; i < n; i++) {
            if(i % 100 == 0)System.out.println("determinantModFast loop i=" + i + " by " + n);
            // 寻找主元
            int piv = i;
            while (piv < n && m[piv * n + i] == 0) {
                piv++;
            }
            if (piv == n) return BigInteger.ZERO; // 列全为0，奇异矩阵

            // 交换行
            if (piv != i) {
                for (int x = 0; x < n; x++) {
                    long tmp = m[i * n + x];
                    m[i * n + x] = m[piv * n + x];
                    m[piv * n + x] = tmp;
                }
                sign = -sign;
            }

            long pivot = m[i * n + i];

            // det = det * pivot mod q
            det = (det * pivot) % q;

            // 计算逆元 (这里用 BigInteger 算一次逆元开销很小)
            // 也可以手写 exgcd 纯 long 实现，但 BigInteger.modInverse 够快了
            long inv = BigInteger.valueOf(pivot).modInverse(BigInteger.valueOf(q)).longValue();

            // 归一化当前行
            for (int j = i; j < n; j++) {
                m[i * n + j] = (m[i * n + j] * inv) % q;
            }

            // 消元
            for (int r = i + 1; r < n; r++) {
                long f = m[r * n + i];
                if (f == 0) continue;

                for (int c = i; c < n; c++) {
                    // m[r][c] = m[r][c] - f * m[i][c]
                    long sub = (f * m[i * n + c]) % q;
                    m[r * n + c] = m[r * n + c] - sub;
                    if (m[r * n + c] < 0) m[r * n + c] += q;
                }
            }
        }

        if (sign == -1) {
            det = -det;
            if (det < 0) det += q;
        }

        return BigInteger.valueOf(det);
    }

    // === 通用通道 (原本的 BigInteger 实现) ===
    // 稍微优化了 GCD 检查，假设 q 为素数可跳过
    private static BigInteger determinantModGeneral(BigIntMatrix1 A) {
        int n = A.getRowDimension();
        if (n != A.getColumnDimension()) return BigInteger.ZERO;
        BigInteger q = A.getModulus();

        BigInteger[] m = new BigInteger[n * n];
        for (int i = 0; i < n * n; i++) {
            // 扁平化访问 A.data (假设你能访问到 data，否则用 get)
            int r = i / n;
            int c = i % n;
            m[i] = A.get(r, c).mod(q);
        }

        BigInteger det = BigInteger.ONE;
        int sign = 1;
        // 简单判断素数优化
        boolean isPrime = q.isProbablePrime(20);

        for (int i = 0; i < n; i++) {
            int piv = i;
            while (piv < n) {
                BigInteger val = m[piv * n + i];
                // 如果是素数，只要非0即互质；否则才算 gcd
                if (!val.equals(BigInteger.ZERO)) {
                    if (isPrime || val.gcd(q).equals(BigInteger.ONE)) break;
                }
                piv++;
            }
            if (piv == n) return BigInteger.ZERO;

            if (piv != i) {
                for (int x = 0; x < n; x++) {
                    BigInteger tmp = m[i * n + x];
                    m[i * n + x] = m[piv * n + x];
                    m[piv * n + x] = tmp;
                }
                sign = -sign;
            }

            BigInteger pivot = m[i * n + i];
            det = det.multiply(pivot).mod(q);
            BigInteger inv = pivot.modInverse(q);

            // 优化：提前判断 inv 是否为 1
            if (!inv.equals(BigInteger.ONE)) {
                for (int j = i; j < n; j++) m[i * n + j] = m[i * n + j].multiply(inv).mod(q);
            }

            for (int r = i + 1; r < n; r++) {
                BigInteger f = m[r * n + i];
                if (f.equals(BigInteger.ZERO)) continue;
                for (int c = i; c < n; c++) {
                    m[r * n + c] = m[r * n + c].subtract(f.multiply(m[i * n + c])).mod(q);
                    // 只有结果为负时才加 q，避免 mod 运算
                    if (m[r * n + c].signum() < 0) m[r * n + c] = m[r * n + c].add(q);
                }
            }
        }
        if (sign == -1) det = det.negate().mod(q);
        return det;
    }

    public static BigInteger GSN_Length(BigIntMatrix1 matrix) {
        BigInteger maxLen = BigInteger.ZERO;
        for (int j = 0; j < matrix.cols; j++) {
            BigInteger sumSq = BigInteger.ZERO;
            for (int i = 0; i < matrix.rows; i++) {
                BigInteger val = matrix.get(i, j);
                sumSq = sumSq.add(val.multiply(val));
            }
            BigInteger colNorm = sqrt(sumSq);
            if (colNorm.compareTo(maxLen) > 0) maxLen = colNorm;
        }
        return maxLen;
    }

    public static BigInteger[] multiplyMatrixVector(BigIntMatrix1 M, BigInteger[] v) {
        int rows = M.rows;
        int cols = M.cols;
        BigInteger q = M.modulus;
        BigInteger[] result = new BigInteger[rows];

        IntStream.range(0, rows).parallel().forEach(i -> {
            BigInteger sum = BigInteger.ZERO;
            for (int j = 0; j < cols; j++) {
                sum = sum.add(M.get(i, j).multiply(v[j]));
            }
            result[i] = sum.mod(q);
        });
        return result;
    }

    // --- Helpers ---

    public void swapCols(int i, int j) {
        for (int k = 0; k < rows; k++) {
            BigInteger temp = data[k * cols + i];
            data[k * cols + i] = data[k * cols + j];
            data[k * cols + j] = temp;
        }
    }

    public void colSub(int i, int j, BigInteger factor) {
        for (int k = 0; k < rows; k++) {
            BigInteger val = data[k * cols + j].multiply(factor);
            data[k * cols + i] = data[k * cols + i].subtract(val);
        }
    }

    public BigIntMatrix1 copy() {
        BigIntMatrix1 c = new BigIntMatrix1(rows, cols, modulus);
        System.arraycopy(data, 0, c.data, 0, data.length);
        return c;
    }

    public static BigIntMatrix1 rowConcat(BigIntMatrix1 A, BigIntMatrix1 B) {
        if (A.cols != B.cols) throw new IllegalArgumentException("Cols mismatch");
        BigIntMatrix1 res = new BigIntMatrix1(A.rows + B.rows, A.cols, A.modulus);
        System.arraycopy(A.data, 0, res.data, 0, A.data.length);
        System.arraycopy(B.data, 0, res.data, A.data.length, B.data.length);
        return res;
    }

    public static BigIntMatrix1 columnConcat(BigIntMatrix1 A, BigIntMatrix1 B) {
        if (A.rows != B.rows) throw new IllegalArgumentException("Rows mismatch");
        BigIntMatrix1 res = new BigIntMatrix1(A.rows, A.cols + B.cols, A.modulus);
        IntStream.range(0, A.rows).parallel().forEach(i -> {
            System.arraycopy(A.data, i * A.cols, res.data, i * res.cols, A.cols);
            System.arraycopy(B.data, i * B.cols, res.data, i * res.cols + A.cols, B.cols);
        });
        return res;
    }

    private void checkSameShape(BigIntMatrix1 other) {
        if (rows != other.rows || cols != other.cols) throw new IllegalArgumentException("Shape mismatch");
    }

    private void swapRows(int i, int j) {
        BigInteger[] tmp = new BigInteger[cols];
        System.arraycopy(data, i * cols, tmp, 0, cols);
        System.arraycopy(data, j * cols, data, i * cols, cols);
        System.arraycopy(tmp, 0, data, j * cols, cols);
    }

    private static BigInteger sqrt(BigInteger n) {
        if (n.signum() < 0) throw new IllegalArgumentException("Negative argument");
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

    // Internal helper for fast reduction
    private static BigInteger fastReduce192ToBigInt(long hi, long mid, long lo, BigInteger modulus) {
        byte[] bytes = new byte[24];
        longToBytes(hi, bytes, 0);
        longToBytes(mid, bytes, 8);
        longToBytes(lo, bytes, 16);
        return new BigInteger(bytes).mod(modulus);
    }

    private static void longToBytes(long val, byte[] buffer, int offset) {
        buffer[offset] = (byte) (val >>> 56);
        buffer[offset + 1] = (byte) (val >>> 48);
        buffer[offset + 2] = (byte) (val >>> 40);
        buffer[offset + 3] = (byte) (val >>> 32);
        buffer[offset + 4] = (byte) (val >>> 24);
        buffer[offset + 5] = (byte) (val >>> 16);
        buffer[offset + 6] = (byte) (val >>> 8);
        buffer[offset + 7] = (byte) (val);
    }

    // --- Legacy Compatibility Methods ---

    @Override
    public String toString() {
        // Warning: Dumping large matrix is slow
        StringBuilder sb = new StringBuilder();
        if (rows > 20 || cols > 20) return "BigIntMatrix[" + rows + "x" + cols + "] (Content hidden)";
        for (int i = 0; i < rows; i++) {
            sb.append("[");
            for(int j=0; j<cols; j++) {
                sb.append(data[i*cols+j]).append(j==cols-1?"":", ");
            }
            sb.append("]").append(System.lineSeparator());
        }
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BigIntMatrix1)) return false;
        BigIntMatrix1 that = (BigIntMatrix1) o;
        if (this.rows != that.rows || this.cols != that.cols) return false;
        // Check exact data match
        return Arrays.equals(this.data, that.data);
    }

    // Original returned String, keeping compatibility
    public String isZeroMatrix() {
        for (BigInteger v : data) {
            if (!v.mod(modulus).equals(BigInteger.ZERO)) return "Not a zero matrix.";
        }
        return "Is a zero matrix.";
    }

    // Added useful boolean version
    public boolean isZeroMatrixBoolean() {
        for (BigInteger v : data) {
            if (!v.mod(modulus).equals(BigInteger.ZERO)) return false;
        }
        return true;
    }
}