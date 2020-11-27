/**
 * SIR model.
 */
model sir {
  const deltat = 0.01; // time step
  const epsilon = 2.0e-6; // observation error tolerance

  param beta;
  param gamma;
  state s, i, r, c;  // susceptible, infectious, recovered, cases
  noise n_r, n_c;
  obs y_c;  // observations

  sub parameter {
    beta ~ gamma(2.0, 1.0);
    gamma ~ gamma(2.0, 1.0);
  }

  sub proposal_parameter {
    beta ~ truncated_normal(beta, 0.01, 0.0);
    gamma ~ truncated_normal(beta, 0.01, 0.0);
  }

  sub initial {
    s <- 1000000;
    i <- 10;
    r <- 0;
    c <- 0;
    n_r <- 0;
  }

  sub transition() {
    inline n = s + i + r;
    n_c ~ binomial(s, 1.0-exp(-beta * i / n * deltat));
    n_r ~ binomial(i, 1.0-exp(-gamma * deltat));
    s <- s - c;
    i <- i + c - n_r;
    r <- r + n_r;
    c <- n_c;
  }

  sub observation {
    y_c ~ poisson(c + epsilon);
  }
}

