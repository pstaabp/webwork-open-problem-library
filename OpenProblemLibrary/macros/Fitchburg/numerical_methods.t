

use strict;
use warnings;


my $webwork_dir = "";
my $pg_dir = "";
my $library_dir = "";

BEGIN {
  $ENV{MOD_PERL_API_VERSION}=2;  # ensure that mod_perl2 is used.
  $webwork_dir = $ENV{WEBWORK_ROOT} || die "The environment variable WEBWORK_ROOT needs to be defined.";
  $pg_dir = $ENV{PG_ROOT};

  $library_dir = "$webwork_dir/../libraries";
  if (not defined $pg_dir) {
    $pg_dir = "$webwork_dir/../pg";
  }

  die "The directory $webwork_dir does not exist" if (not -d $webwork_dir);
  die "The directory $pg_dir does not exist" if (not -d $pg_dir);

}

use lib "$webwork_dir/lib";
use lib "$pg_dir/lib";
use lib "$pg_dir/macros";
use lib ".";

use Test2::Bundle::More;
use Test2::Tools::Compare qw/is/;
#use Test2::Exception;

use Parser;
use Value;
use Class::Accessor;
use PGcore;
use Try::Tiny;
use Data::Dump qw/dd/;

$WeBWorK::run_test = 1;
#
#  Fake these for now
#
sub loadMacros {
  foreach my $file (@_) {
    require $file;
    my $init = $file;
    $init =~ s!.*/!!;
    $init =~ s!\.[^.]*$!!;
    $init = "_${init}_init";
    eval("$init()");
  }
}

sub be_strict { use strict;}

sub ParserDefineLog {eval {sub log($) {CommonFunction->Call("log",@_)}}};

require "Value.pl";
require "Parser.pl";

my %context = ();

sub Context {Parser::Context->current(\%context,@_)}
unless (%context && $context{current}) {
  # ^variable our %context
  %context = ();  # Locally defined contexts, including 'current' context
  # ^uses Context
  Context();      # Initialize context (for persistent mod_perl)
}

require 'numericalMethods.pl';


Context("Numeric");
# Context()->flags->set(tolerance => 1e-10);

### testing the round function
subtest 'Rounding tests' => sub {
  is(round_digits(1.23456,3),1.23,'rounding to 3 digits');
  is(round_digits(1.23456,5),1.2346,'rounding to 5 digits');
  is(round_digits(123.456,4),123.5,'rounding to 4 digits');
};


my $f = Compute("x^2-2");

my @known;
$known[0][0]=0; $known[0][1] =2;
$known[1][0]=1; $known[1][1] =2;
$known[2][0]=1; $known[2][1] =1.5;
$known[3][0]=1.25; $known[3][1] =1.5;
$known[4][0]=1.375; $known[4][1] =1.5;
$known[5][0]=1.375; $known[5][1] =1.4375;
$known[6][0]=1.40625; $known[6][1] =1.4375;
$known[7][0]=1.40625; $known[7][1] =1.421875;
$known[8][0]=1.4140625; $known[8][1] =1.421875;
$known[9][0]=1.4140625; $known[9][1] =1.41796875;
#$known[9][0]=1.4140625; $known[9][1] =1.416015625;

my @ints = bisection($f,0.0,2.0,1e-12,9);

subtest 'Bisection tests' => sub {
  is(\@ints,\@known);

  ## create a function without a root in the given interval

  my $err;
  my $g = Compute("x^2+2");
  try {
    my @out = bisection($g,0.0,2.0);
  } catch {
    $err = $_;
    #warn "caught error: $_"; # not $@
  };
  is(substr($err,0,20),substr('To ensure there is a root, f(a)f(b) must be less than 0',0,20), "this function doesn't have a root");
};

### Newton's method

my @new_seq = (20.,10.05000000,5.124502488,2.757392139,1.741357581,1.444943382,1.414540330);

my @x = newtonsMethod($f,20.0,1e-12,5);

subtest "Newton's Method tests" => sub {
  is(\@x,\@new_seq);
};

subtest "Newton Divided Difference Test" => sub {

  my $x = [0, 0.785398, 1.5708, 2.35619, 3.14159];
  my $y = [0, 1, 0, -1, 0];

  warn join(", ", @$x);
  warn join(", ", @$y);

  my @a = newton_divided_difference($x,$y);
};


done_testing();
