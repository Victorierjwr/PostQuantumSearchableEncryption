package qut.edu.au.pqcsebase.algorithms;

import java.math.BigInteger;
import java.security.SecureRandom;

public class SampleZ {

    private static SecureRandom random = new SecureRandom();

    /**
     * Sample z from discrete Gaussian distribution D_{Z, sigma, c}
     *
     * @param sigma standard deviation sigma = n^c
     * @param c center c > 0
     * @param n dimension
     * @return sampled integer x
     */
    public static BigInteger sampleZ(double sigma, double c, long n) {

//        double log2n = Math.log(n) / Math.log(2);
//        double t = Math.sqrt(log2n);

        // t(n) = log2(n)
        double t = Math.log(n) / Math.log(2);

        double z = sigma * t;
        double f = c - z; // 下界 (float)
        double b = c + z; // 上界 (float)

        long min = (long) Math.ceil(f);
        long max = (long) Math.floor(b);

        long rangeSize = max - min + 1;
        if (rangeSize <= 0) {
            return BigInteger.valueOf(Math.round(c));
        }

        while(true){
            long y = min + (long)(random.nextDouble() * rangeSize);

            double dist = y - c;               // (x - c)
            double top = Math.PI * dist * dist; // π * (x - c)^2
            double bottom = sigma * sigma;             // sigma^2
            double g = top / bottom;           // π(x-c)^2 / sigma^2
            double prob = Math.exp(-g);        // exp( - π(x-c)^2 / sigma^2 )

            if (random.nextDouble() < prob) {
                return BigInteger.valueOf(y);
            }
        }
    }
}
