package qut.edu.au.pqcsebaseopt.projects.IBEKS.vo;

import lombok.*;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

/**
 * Master Secret Key class for IBEKS scheme.
 */
@Data@Setter@Getter@AllArgsConstructor@ToString
public class MSK {

    // m*m dimension matrix basis
    public BigIntMatrix msk;
}
