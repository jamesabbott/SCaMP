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
use File::Basename;
use File::Copy qw(copy);
use File::Copy::Recursive qw(dircopy);
use File::Path qw(make_path);

use vars '$AUTOLOAD';

use warnings;
use strict;

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

    my ( $self, %args ) = @_;
    my $scamp_root = $args{'scamp_root'} || croak "scamp_root argument is not prorvided";
    $self->{"_scamp_root"} = $scamp_root;

    my $config = LoadFile("$scamp_root/etc/SCaMP.yaml");

    foreach my $key ( keys %$config ) {
        $self->{"_${key}"} = $config->{$key};
    }

    return ();
}

=pod
 
=over
 
=item B<get_job_id>
 
    Returns the id of the task being executed. This is derived from the setting
    of the PBS_JOBID or JOB_ID  environmental variables which are set
    during the execution of an array job task. 

=back
 
=cut

sub get_job_id {

    my $self = shift;
    my $job;

    if ( exists( $ENV{'JOB_ID'} ) ) {
        $job = $ENV{'JOB_ID'};
    }
    elsif ( exists( $ENV{'PBS_JOBID'} ) ) {
        $job = $ENV{'PBS_JOBID'};
    }
    else {
        $job = 1;
    }
    return ($job);
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

    if ( exists( $ENV{'SGE_TASK_ID'} ) ) {
        $task = $ENV{'SGE_TASK_ID'};
    }
    elsif ( exists( $ENV{'PBS_ARRAY_INDEX'} ) ) {
        $task = $ENV{'PBS_ARRAY_INDEX'};
    }
    else {
        $task = 1;
    }
    return ($task);
}

=pod

=over 

=item B<setup_paths>

  Create necessary directories for run, staging data to TMP_DIR if requested. 
  Required arguments
    $ ('stage' - boolean: should data be staged to local storage?)
    $ ('sample' - sample name)
    $ ('src_dir' - origin dir of data)
    $ ('work_dir' - location for temp working space)
    $ ('job_files' - arrayref of filename to stage)

  Returns:
    $ (in_dir)
    $ (scratch_dir)

=back

=cut

sub setup_paths {

    my ( $self, %args ) = @_;
    my $stage = $args{'stage'};                                          ## can be passed 0 which fails assignment check
    my $sample = $args{'sample'} || croak "stage argument not provided";
    my $src_dir   = $args{'src_dir'} || croak "src_dir argument not provided";
    my $work_dir  = $args{'work_dir'} || croak "work_dir argument not provided";
    my $job_files = $args{'job_files'} || croak "job_files argument not provided";
    my $db = $args{'db'}; 
    my $db_dir = $self->get_database_dir();

    my ( $in_dir, $scratch_dir );

    # for staged data, create local tmp dir in $TMP_DIR, copy input files to $TMP_DIR
    if ($stage) {
        $in_dir      = $ENV{'TMPDIR'};
        $scratch_dir = "$ENV{'TMPDIR'}/work";
        mkdir "$scratch_dir" or croak "Error creating $scratch_dir: $!";
        mkdir "$in_dir/db" or croak "Error creating $in_dir/db: $!";

        print "Staging data to $in_dir...\n";
        foreach my $file (@$job_files) {
            my $basename = fileparse($file);
            print "copy $file ->  $in_dir/$basename\n";
            copy( $file, "$in_dir/$basename" ) or die "Error copying $file -> $in_dir/$basename:$!";
        }

	if ($db) {
	    print "copying $db database -> $in_dir/db\n";
	    my $resolved=readlink("$db_dir/$db/latest") or die "Error reading $db_dir/$db/lastest symlink: $!";
	    dircopy("$db_dir/$db/$resolved", "$in_dir/db") or croak "Error copying $db: $!";
	    $db_dir = "$in_dir/db";
	}
    }
    else {
        $in_dir      = $src_dir;
        $scratch_dir = "$work_dir/tmp/$sample";
        if ( !-d $scratch_dir ) {
            make_path( "$scratch_dir", { error => \my $err } );
            if (@$err) {
                for my $diag (@$err) {
                    my ( $file, $message ) = %$diag;
                    if ( $file eq '' ) {
                        print "general error: $message\n";
                    }
                    else {
                        print "problem unlinking $file: $message\n";
                    }
                }
            }
        }
    }
    print "in_dir = $in_dir, scratch_dir = $scratch_dir, db_dir = $db_dir\n";
    
    return ( $in_dir, $scratch_dir, $db_dir );
}

1;
