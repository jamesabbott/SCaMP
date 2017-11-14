package SCaMP;

=pod

=head1 NAME

    SCaMP - provides centralised config data

=head1 SYNOPSIS

    my $scamp = SCaMP->new();
    my $database_dir = $scamp->get_database_dir();

=head1 DESCRIPTION

    Provides common methods and configuration data for SCaMP.
    Configuration data is loaded from etc/SCaMP.yaml.
    An autoloaded accessor is created for each top-level entry in the YAML file 
    i.e. get_databases()

=cut

use Carp qw(croak);
use YAML::XS qw(LoadFile);

use vars '$AUTOLOAD';

{
    my %_attrs = (
                   _database_dir => ['read'],
                   _databases    => ['read'],
                   _work_dir     => ['read'],
                   _scratch_dir  => ['read'],
                 );

    sub _accessible {
        my ( $self, $attr, $mode ) = @_;
        return $_attrs{$attr}[0] =~ /$mode/;
    }

}

sub DESTROY {
    my $self = shift;
    return ();
}

sub AUTOLOAD {

    my ( $self, $newval ) = @_;

    $AUTOLOAD =~ /.*::get(_\w+)/
      and $self->_accessible( $1, 'read' )
      and return $self->{$1};

    croak "No such method :$AUTOLOAD";

}

sub new {

    my ( $class, %args ) = @_;
    my $self = {};
    bless $self, $class;
    $self->_init(%args);

    return $self;

}

sub _init {

    my $self   = shift;
    my $scamp_root = shift;
    $self->{"_scamp_root"} = $scamp_root;

#    my $config = LoadFile("$FindBin::Bin/../etc/SCaMP.yaml");
    my $config = LoadFile("$scamp_root/etc/SCaMP.yaml");

    foreach my $key ( keys %$config ) {
        $self->{"_${key}"} = $config->{$key};
    }

    return ();
}

=pod
 
=over
 
=item B<get_task_id>
 
    Returns the id of the task being executed. This is derived from the setting
    of the SGE_TASK_ID or PBS_ARRAY_INDEX environmental variables which are set
    during the execution of an array job task. 

=back
 
=cut

sub get_task_id {

    my $self = shift;
    my $task;

    if ( exists($ENV{'SGE_TASK_ID'}) ) {
        $task = $ENV{'SGE_TASK_ID'};
    }
    elsif ( exists($ENV{'PBS_ARRAY_INDEX'}) ) {
        $task = $ENV{'PBS_ARRAY_INDEX'};
    }
    else {
        croak "This script should be run as an array job using either the SGE or PBSPro batch queueing systems";
    }
    return ($task);
}

1;
