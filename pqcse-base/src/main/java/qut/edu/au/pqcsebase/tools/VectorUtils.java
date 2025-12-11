package qut.edu.au.pqcsebase.tools;


import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Objects;

public class VectorUtils {

    private static final SecureRandom RANDOM = new SecureRandom();

    /**
     * 生成随机向量，元素均在 [0, q-1]
     */
    public static BigInteger[] randomVector(int length, BigInteger q) {
        if (length <= 0) throw new IllegalArgumentException("length must be > 0");
        Objects.requireNonNull(q, "q must not be null");
        if (q.compareTo(BigInteger.ONE) <= 0) throw new IllegalArgumentException("q must be > 1");

        BigInteger[] v = new BigInteger[length];
        for (int i = 0; i < length; i++) {
            BigInteger val;
            // 生成 [0, q-1]，避免偏差
            do {
                val = new BigInteger(q.bitLength(), RANDOM);
            } while (val.compareTo(q) >= 0);
            v[i] = val;
        }
        return v;
    }

    /**
     * 普通整数内积 (不取模)
     */
    public static BigInteger innerProduct(BigInteger[] a, BigInteger[] b) {
        validateSameLength(a, b);
        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < a.length; i++) {
            sum = sum.add(a[i].multiply(b[i]));
        }
        return sum;
    }

    /**
     * 模 q 的内积，结果在 [0, q-1]
     */
    public static BigInteger innerProductModQ(BigInteger[] a, BigInteger[] b, BigInteger q) {
        validateSameLength(a, b);
        Objects.requireNonNull(q, "q must not be null");
        if (q.compareTo(BigInteger.ONE) <= 0) throw new IllegalArgumentException("q must be > 1");

        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < a.length; i++) {
            BigInteger ai = a[i].mod(q);
            BigInteger bi = b[i].mod(q);
            sum = sum.add(ai.multiply(bi)).mod(q);
        }
        // 保证非负
        if (sum.signum() < 0) sum = sum.add(q);
        return sum;
    }

    private static void validateSameLength(BigInteger[] a, BigInteger[] b) {
        Objects.requireNonNull(a, "a must not be null");
        Objects.requireNonNull(b, "b must not be null");
        if (a.length != b.length) throw new IllegalArgumentException("Vectors must have same length");
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null || b[i] == null) throw new IllegalArgumentException("Vector entries must not be null");
        }
    }
}

