model corrMod {
  param sigPN;
  param sigOE;
  param a;
  param b;
  noise w;
  state x;
  obs y;

  sub parameter {
    sigPN ~ uniform(0.0001,1);
    sigOE ~ uniform(0.0001,1);
    a ~ uniform(-.9999,.9999);
    b ~ gaussian(0, 1000);
  }

sub lookahead_transition {
   x <-  0.95*x + 1;
}

  sub initial {
    x ~ gaussian(mean = b/(1-a),std = sqrt(pow(sigPN,2)/(1-pow(a,2))));
  }

  sub transition {
    w ~ gaussian(mean=0, std=sigPN);
    x <- a*x + b + w;
  }
	
  sub observation {
    y ~ gaussian(mean=x, std=sigOE);
  }
}
