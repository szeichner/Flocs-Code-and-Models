function[m1, m2] = CalculateLeastSquares(x, y)
    G = [ones(length(x), 1) x(:)];
    m_lls = (G'*G)^(-1) * G' * y;
    m1 = m_lls(1);
    m2 = m_lls(2);  
end