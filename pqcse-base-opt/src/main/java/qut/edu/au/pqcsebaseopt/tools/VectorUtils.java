package qut.edu.au.pqcsebaseopt.tools;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Objects;
import java.util.concurrent.ThreadLocalRandom;

public class VectorUtils {

    private static final ThreadLocalRandom RANDOM = ThreadLocalRandom.current();

    /**
     * 生成随机向量，元素均在 [0, q-1]
     */
    public static BigInteger[] randomVector(int length, BigInteger q) {
        if (length <= 0) throw new IllegalArgumentException("length must be > 0");
        Objects.requireNonNull(q, "q must not be null");
        if (q.compareTo(BigInteger.ONE) <= 0) throw new IllegalArgumentException("q must be > 1");

        BigInteger[] v = new BigInteger[length];
        int numBits = q.bitLength();

        for (int i = 0; i < length; i++) {
            BigInteger val;
            // 生成 [0, q-1]，避免偏差
            do {
                val = new BigInteger(numBits, RANDOM);
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
     * 优化版：模 q 的内积，结果在 [0, q-1]
     * 优化点：移除循环内的 mod 操作，只在最后取模。
     */
    public static BigInteger innerProductModQ(BigInteger[] a, BigInteger[] b, BigInteger q) {
        validateSameLength(a, b);
        Objects.requireNonNull(q, "q must not be null");
        if (q.compareTo(BigInteger.ONE) <= 0) throw new IllegalArgumentException("q must be > 1");

        BigInteger sum = BigInteger.ZERO;

        // 核心优化：累加过程中不取模
        // BigInteger 可以处理任意大整数，加法和乘法比取模(除法)快得多。
        for (int i = 0; i < a.length; i++) {
            sum = sum.add(a[i].multiply(b[i]));
        }

        // 最后统一取模。BigInteger.mod 保证结果非负。
        return sum.mod(q);
    }

    public static BigInteger[] addVectorsMod(BigInteger[] vector1, BigInteger[] vector2, BigInteger q) {
        validateSameLength(vector1, vector2); // 添加长度检查更安全
        BigInteger[] result = new BigInteger[vector1.length];
        for (int i = 0; i < vector1.length; i++) {
            // 这里每次必须取模，防止结果无限增长，且通常作为中间结果需要保持在域内
            result[i] = vector1[i].add(vector2[i]).mod(q);
        }
        return result;
    }

    /**
     * 优化版验证：只检查引用和长度，移除 O(N) 的 null 检查
     */
    private static void validateSameLength(BigInteger[] a, BigInteger[] b) {
        if (a == null || b == null) {
            throw new NullPointerException("Vectors must not be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("Vectors must have same length: " + a.length + " != " + b.length);
        }
    }
}