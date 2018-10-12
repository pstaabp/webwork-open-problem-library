#use feature "say";

#use Data::Dump qw/dump dd/;

# These macros handle many basic numerical methods

# the following is needed if the test suite is run on this file.

if ($WeBWorK::run_test){
  warn $WeBWorK::run_test;
  loadMacros(
    "MathObjects.pl",
    "PGauxiliaryFunctions.pl"
  );

  #use Data::Dump qw/dd/;
}



####
#
# the subroutine newtonsMethod applies newton's method for the given function
# the method will stop when |x_{n+1}-x_n| < $eps or if it reaches $max steps
#
# inputs:
#   $f -- a Mathobject as a function of x
#   $x0 -- the initial point (a number)
#   $eps -- the stopping criteria.  If this isn't defined, then use 1e-6
#   $max -- the maximum number of steps to use.  if this isn't defined, use 20
#
# output:
#   an array of x-values.
#
#######

sub newtonsMethod {
  my ($f,$x0,$eps,$max) = @_;


  $f->perlFunction("F");
  $df = $f->D('x');
  $df->perlFunction("dF");
  $eps = 1e-6 unless defined($eps);
  $max = 20 unless defined($max);

  @x=($x0);
  $n = 0;
  LOOP: {
    do {
      my $f0 =F($x[$n]);
      my $df0 = dF($x[$n]);
      my $x1 = $x[$n]-$f0/$df0;
      push(@x,$x[$n]-$f->eval(x=>$x[$n])/$df->eval(x=>$x[$n]));
      if ($n>=$max){ last LOOP;}
      $n++;
    } while (abs($x[$n]-$x[$n-1])>$eps);
  }
  return @x;
}

####
#
#  the subroutine bisection performs the bisection method on a input function
#  and left and right endpoints.
#
#  inputs:
#     f: a MathObject function
#     a: the initial left endpoint
#     b: the initial right endpoint
#     $eps -- the stopping criteria.  If this isn't defined, then use 1e-6
#     $max -- the maximum number of steps to use.  if this isn't defined, use 20
#
#  output:
#     an array of intervals (each interval is a length-2 array)
#
####

sub bisection {
  my ($f,$a0,$b0,$eps,$max) = @_;
  my $fa = $f->eval(x => $a0);
  my $fb = $f->eval(x => $b0);
  if ($fa*$fb>0 ) {
    die "To ensure there is a root, f(a)f(b) must be less than 0";
  }

  $eps = 1e-6 unless defined($eps);
  $max = 20 unless defined($max);
  my @ints;
  $ints[0][0] = $a0; $ints[0][1]=$b0;

  for($i=0;$i<$max;$i++){
    my $c = 0.5*($ints[$i][0]+$ints[$i][1]);
    my $fa = $f->eval(x => $ints[$i][0]);
    my $fb = $f->eval(x => $ints[$i][1]);
    my $fc = $f->eval(x => $c);
    $ints[$i+1][0] =($fa*$fc<0)?$ints[$i][0]:$c;
    $ints[$i+1][1] =($fa*$fc<0)?$c:$ints[$i][1];
  }
  return @ints;
}

####
#
# the subroutine secant performs the secant method on the given function to the
# given number of digits.
#
#  inputs:
#     f: a MathObject function
#     x0: the first initial point
#     x1: the second initial point
#     $eps -- the stopping criteria.  If this isn't defined, then use 1e-6
#     $max -- the maximum number of steps to use.  if this isn't defined, use 20
#
#  output:
#     an array of intervals (each interval is a length-2 array)
#
####

sub secant {
  my ($f,$x0,$x1, $eps,$max) = @_;


  $f->perlFunction("F");
  $eps = 1e-6 unless defined($eps);
  $max = 20 unless defined($max);

  @x = ($x0,$x1);

  for($i=2;$i<=$max;$i++){
    $m=($f->eval(x=>$x[$i-2])-$f->eval(x=>$x[$i-1]))/($x[$i-2]-$x[$i-1]);
    push(@x,$x[$i-1]-$f->eval(x=>$x[$i-1])/$m);
    if(abs($x[$i]-$x[$i-1])<$eps){
      return @x;
    }
  }

  return @x;
}

#####
#
# the subroutine newton_divided_difference creates a
# Newton Divided Difference Table and returns it as an 2D array
#
#
#  inputs:
#    x: an array of x values
#    y: an array of y values corresponding to the x values.
#
#  output:
#    an array representing the divided difference table.
#
#  note if @a is the array
#
#  x[0]  y[0]
#              a[0][0]
#  x[1]  y[1]           a[1][1]
#              a[0][1]           a[2][2]
#  x[2]  y[2]           a[1][2]
#              a[0][2]
#  x[3]  y[3]
#
#
#####

sub newton_divided_difference {
  my ($x,$y) = @_;

  if (ref($x) eq 'ARRAY'){
     @x = @$x;
  } elsif (ref($$x) eq 'ARRAY'){
    @x = @$$x;
  }
  if (ref($y) eq 'ARRAY'){
    @y = @$y;
  } elsif (ref($$x) eq 'ARRAY'){
    @y = @$$y;
  }

  @x = map {  Value::isValue($_) ? $_->value : $_ } @x;
  @y = map {  Value::isValue($_) ? $_->value : $_ } @y;

  die "The lengths of the input arrays must be the same" unless (scalar(@x) == scalar(@y));

  my $n = scalar(@x);

  my @a = ();

  for(my $j=0;$j<$n-1;$j++){
    @a[$j] = ();
    for(my $i=0; $i<$n-$j-1; $i++){
      $a[$j][$i+$j] = ($j==0)? ($y[$i+1]-$y[$i])/($x[$i+1]-$x[$i])
                        : ($a[$j-1][$i+$j]-$a[$j-1][$i+$j-1])/($x[$i+$j+1]-$x[$i]);
    }
  }

  return @a;
}

####
#
#  the subroutine round_digits returns the given input to the number of digits.
#  Note: the number of digits is the number of significant figures.
#
#  inputs:
#    x: the given number
#    digits: a positive integer which is the number of digits.
#
#   ouput: x rounded to the given number of digits
#
#  example:
#    round_digits(1.234567,3)  returns 1.23
#
######


sub round_digits {
  ($x,$digits) = @_;
  $n = floor(log(abs($x))/log(10));
  return round($x*10**($digits-1-$n))/10**($digits-1-$n);
}



1;
