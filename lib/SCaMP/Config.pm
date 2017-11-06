package SCaMP::Config;

=pod

=head1 NAME

    SCaMP - provides centralised config data

=head1 SYNOPSIS

    my $scamp_config = SCaMP::Config->new();
    my $config = $scamp_config->get_config(); 

=head1 DESCRIPTION

    Provides configuration data from etc/SCaMP.yaml via an OO interface
    An autoloaded accessor is created for each top-level entyr in the YAML file 
    i.e. get_databases()

=cut

use Carp qw(croak);
use YAML::XS qw(LoadFile);

use vars '$AUTOLOAD';                                

{                                                    
    my %_attrs = (                                   
                   _database_dir    => ['read'],     
                   _databases	    => ['read'],     
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

	my $self = shift;
	my $config = LoadFile("$FindBin::Bin/../etc/SCaMP.yaml");
	foreach my $key(keys %$config) {
	    $self->{"_${key}"}=$config->{$key};
	}

	return();
}

1;
