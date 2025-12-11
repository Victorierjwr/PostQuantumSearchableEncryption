package qut.edu.au.pqcsebase.projects.IBEKS.vo;

import lombok.*;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;

/**
 * Public Parameters class for IBEKS scheme.
 */
@Data
@ToString
@AllArgsConstructor@NoArgsConstructor
@Getter@Setter
public class PP {

    // A n*m mod q matrix
    public BigIntMatrix A;
    // B n*m mod q matrix
    public BigIntMatrix B;

    //n-N-1 dimension vectors
    public BigInteger[] b1;
    public BigInteger[] y1;

    // n dimension vectors
    public BigInteger[] u;
}
