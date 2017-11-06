package SCaMP::Config;

=pod

=head1 NAME

    SCaMP - provides centralised config data

=head1 SYNOPSIS

    my $scamp_config = SCaMP::Config->new();
    my $config = $scamp_config->get_config(); 

=head1 DESCRIPTION

    Provides configuration data from etc/SCaMP.yaml via an OO interface

=cut

use Carp qw(croak);
use YAML::XS qw(LoadFile);

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

	$self->{'config'}=$config;
}

sub get_config {
	my $self = shift;
	return($self->{'config'})
}

1;
