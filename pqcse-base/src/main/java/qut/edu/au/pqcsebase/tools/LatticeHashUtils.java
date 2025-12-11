package qut.edu.au.pqcsebase.tools;

import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import static qut.edu.au.pqcsebase.tools.BigIntMatrix.determinantMod;

/**
 * Implements Hash functions H1 and H2 for Lattice-based schemes.
 * H1: {0,1}* -> Z_q^{mxm} (Invertible Matrix)
 * H2: {0,1}* -> Z_q (Scalar)
 */
public class LatticeHashUtils {

    private final BigInteger q;
    private final int m;

    public LatticeHashUtils(BigInteger q, int m) {
        this.q = q;
        this.m = m;
    }

    /**
     * H2: Maps input string to an integer in Z_q
     */
    public BigInteger StringHashToValueModQ(String input) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hashBytes = digest.digest(input.getBytes(StandardCharsets.UTF_8));
            return new BigInteger(1, hashBytes).mod(this.q);
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-256 algorithm not found", e);
        }
    }

    /**
     * H1: Maps input string to an Invertible Matrix in Z_q^{mxm}
     * Returns a BigIntMatrix object.
     */
    public BigIntMatrix StringHashToMatrixModQ(String input) {
        int attemptCounter = 0;

        while (true) {
            // 1. 生成矩阵 (Deterministic based on input + counter)
            String currentSeed = input + "||" + attemptCounter;
            BigIntMatrix matrix = generateDeterministicMatrix(currentSeed);

            // 2. 检查可逆性
            // 直接调用你的 BigIntMatrix 工具类中的 determinantMod 方法
            BigInteger det = determinantMod(matrix);

            // 只要行列式与 q 互质（若 q 为素数，则 det != 0），则矩阵可逆
            if (!det.equals(BigInteger.ZERO) && det.gcd(q).equals(BigInteger.ONE)) {
                return matrix;
            }

            // 3. 拒绝采样：尝试下一个 nonce
            attemptCounter++;

            if (attemptCounter > 1000) {
                throw new RuntimeException("Failed to find invertible matrix after 1000 attempts.");
            }
        }
    }

    /**
     * Helper: Fill a BigIntMatrix deterministically using SHA-256
     */
    private BigIntMatrix generateDeterministicMatrix(String seed) {
        BigIntMatrix matrix = new BigIntMatrix(m, m, q);
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            int elementCount = 0;

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    // 构造独立的输入以保证每个元素不同但确定
                    String inputForCell = seed + "||" + elementCount;
                    byte[] hashBytes = digest.digest(inputForCell.getBytes(StandardCharsets.UTF_8));

                    BigInteger val = new BigInteger(1, hashBytes).mod(this.q);

                    // 使用 BigIntMatrix 的 setter
                    matrix.set(i, j, val);

                    elementCount++;
                }
            }
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
        return matrix;
    }
}
