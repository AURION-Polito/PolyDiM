from Elliptic_PCC_2D.test_definition import ITest, EllipticPolynomialProblem

def create_test(test_id: int) -> ITest:

    if test_id == 1:
        return EllipticPolynomialProblem()
    else:
        raise ValueError("not valid test id")
