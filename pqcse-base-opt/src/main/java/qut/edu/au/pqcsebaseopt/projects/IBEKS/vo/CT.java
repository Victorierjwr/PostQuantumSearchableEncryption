package qut.edu.au.pqcsebaseopt.projects.IBEKS.vo;

import lombok.*;

import java.math.BigInteger;

@Data
@Getter@Setter
@AllArgsConstructor @NoArgsConstructor
@ToString
public class CT {

    private BigInteger c_u;
    private BigInteger[] c_id;
    private BigInteger[] c_w;

}

