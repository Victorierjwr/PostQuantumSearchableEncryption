package qut.edu.au.pqcsebaseopt.algorithms;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.concurrent.ThreadLocalRandom;

public class SampleZ {

    private static final ThreadLocalRandom random = ThreadLocalRandom.current();

    /**
     * Sample z from discrete Gaussian distribution D_{Z, sigma, c}
     * Optimized for performance.
     */
    public static BigInteger sampleZ(double s, double c, long n) {
        // t(n) = log2(n)
        double t = Math.log(n) / Math.log(2);
        double z = s * t;

        long min = (long) Math.ceil(c - z);
        long max = (long) Math.floor(c + z);
        long rangeSize = max - min + 1;

        if (rangeSize <= 0) {
            return BigInteger.valueOf(Math.round(c));
        }

        while(true){
            // 使用 ThreadLocalRandom 替代 static SecureRandom
            // 解决 StackOverflowError 和多线程竞争问题
            double rand1 = ThreadLocalRandom.current().nextDouble();

            long y = min + (long)(rand1 * rangeSize);

            double dist = y - c;
            double g = (Math.PI * dist * dist) / (s * s);

            // 优化: 快速拒绝
            if (g > 50) continue;

            double prob = Math.exp(-g);

            if (ThreadLocalRandom.current().nextDouble() < prob) {
                return BigInteger.valueOf(y);
            }
        }
    }
}
