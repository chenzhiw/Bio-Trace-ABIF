package Bio::Trace::ABIF;


# A note about unpacking binary values:
#
# According to Applied Biosystems's documentation, some
# data in the ABI files is stored as *signed* short/long.
# In principle, these should be unpack()ed with 's' and 'l',
# respectively. Unfortunately, this seems to work
# only in big-endian architectures. To overcome the problem,
# we use 'n' and 'N' instead, which seems should be portable,
# at least as long as we parse non-negative numbers.
#
# For floating point numbers the situation is worse, because
# there is no portable way of representing floating point numbers
# across different architectures.
#
# Although not documented, float numbers in ABI files
# apparently use standard IEEE representation. So, we adopt the
# following approach: they are decoded as 32 bit strings
# (with descending bit order inside each byte)
# and then they are converted into decimal numbers.
#
use warnings;
use strict;

=head1 NAME

Bio::Trace::ABIF - Perl extension for reading and parsing ABIF (Applied Biosystems, Inc. Format) files
       
=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

The ABIF file format is a binary format for storing data, developed
by Applied Biosystems, Inc. This module provides general methods for accessing
any chunk of information contained into an ABIF file. Besides, it provides
shortcut methods for reading the most commonly used parts of the file, and
methods for manipulating the data (e.g., for computing LOR scores).

  use Bio::Trace::ABIF;
  
  my $ab1 = Bio::Trace::ABIF->new();
  $ab1->open_abif('/Path/to/my/file.ab1');
  
  print $ab1->sample_name(), "\n";
  my @quality_values = $ab1->quality_values();
  my $sequence = $ab1->sequence();
  # etc...

  $ab1->close_abif();
  
If you cannot find a method to retrieve the information you need, you may
use get_data_item() or get_directory().

=cut

#use 5.008006;
use Carp;

# Another useless comment

require Exporter;

our @ISA = qw(Exporter);
my $Debugging = 0;
my $DIR_SIZE = 28; # Size, in bytes, of a directory entry in an ABIF file

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Bio::Trace::ABIF ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

# Standard types
our %TYPES = (
		1 => 'byte', 2 => 'char', 3 => 'word', 4 => 'short', 5 => 'long',
		7 => 'float', 8 => 'double', 10 => 'date', 11 => 'time', 18 => 'pString',
		19 => 'cString', 12 => 'thumb', 13 => 'bool', 6 => 'rational', 9 => 'BCD',
		14 => 'point', 15 => 'rect', 16 => 'vPoint', 17 => 'vRect', 20 => 'tag',
		128 => 'deltaComp', 256 => 'LZWComp', 384 => 'deltaLZW'
	); # User defined data types have tags numbers >= 1024


=head1 CONSTRUCTOR

=head2 new()

  Usage    : my $ab1 = Bio::Trace::ABIF->new();
  Returns  : an instance of ABIF.

Creates an ABIF object.

=cut
	
sub new {
	my $class = shift;
	my $foo = shift;
	my $self = {};
	$self->{'_FH'} = undef; # ABIF file handle
	$self->{'_NUMELEM'} = undef;
	$self->{'_DATAOFFSET'} = undef;
	# Data type codes as specified by AB specification
	$self->{'TYPES'} = \%TYPES;
	bless($self, $class);
	
	return $self;
}

=head1 METHODS

=head2 open_abif()

  Usage    : $ab1->open_abif($pathname) or die 'Wrong filename or format';
  Returns  : 1 if the file is successfully opened and it is an ABIF file;
             0 otherwise. 

Opens the specified file in binary format and checks whether it is in ABIF format.         

=cut

sub open_abif {
	my $self = shift;
	my $filename = shift;
	# Close previously opened file, if any
	if (defined($self->{'_FH'})) { $self->close_abif(); }
	open($self->{'_FH'}, "<", $filename) or return 0;
	binmode($self->{'_FH'});
	unless ($self->is_abif_format()) {
		print STDERR "$filename is not an AB file...\n";
		close($self->{'_FH'});
		$self->{'_FH'} = undef;
		return 0;
	}
	# Determine the number of items (stored in bytes 18-21)
 	# and the offset of the data (stored in bytes 26-29)
	my $bytes;
	unless (seek($self->{'_FH'}, 18, 0)) {
		carp "Error on seeking file $filename";
		return 0;
	}
	# Read bytes 18-29
	unless (read($self->{'_FH'}, $bytes, 12)) {
		carp  "Error on reading $filename";
		return 0;
	}
	# Unpack a 32 bit integer, skip four bytes and unpack another 32 bit integer
	($self->{'_NUMELEM'}, $self->{'_DATAOFFSET'}) = unpack('Nx4N', $bytes);
	
	# Cache tags positions
	$self->_scan_tags();
	
	return 1;
}

# Performs a linear scan of the file,
# and stores the tags's offsets in a hash, for fast retrieval
sub _scan_tags {
	my $self = shift;
	my ($tag_name, $tag_number, $field);
	my $i = 0;
	$self->{'_TAG_INDEX'} = { };
	do {
		my $offset = $self->{'_DATAOFFSET'} + ($DIR_SIZE * $i);
		seek($self->{'_FH'}, $offset, 0);
		# Read Tag name and number (8 bytes)
		read($self->{'_FH'}, $field, 8);
		#($t1, $t2) = unpack('A4l', $field);
		($tag_name, $tag_number) = unpack('A4N', $field);
		${$self->{'_TAG_INDEX'}}{$tag_name . $tag_number} = $offset;
		$i++;
	}
	while ($i < $self->{'_NUMELEM'});  

	return;
}

=head2 close_abif()

  Usage    : $ab1->close_abif();
  Returns  : Nothing.

Closes the currently opened file.

=cut

sub close_abif {
	my $self = shift;
	close($self->{'_FH'});
	foreach my $k (keys %$self) {
		$self->{$k} = undef if $k =~ /^\_/;
	}
}

=head2 is_abif_open()

  Usage    : if ($ab1->is_abif_open()) { # ...
  Returns  : 1 if an ABIF file is open; 0 otherwise.

=cut

sub is_abif_open {
	my $self = shift;
	return defined $self->{'_FH'};
}


=head2 is_abif_format()

  Usage    : if ($ab1->is_abif_format()) { # Ok, it is ABIF
  Returns  : 1 if the file is in ABIF format, 0 otherwise.

Checks that the file is in ABIF format. This method is called automatically
when a file is opened.

=cut

sub is_abif_format {
	my $self = shift;
	my $file_signature;
	# Move to the beginning of the file
	unless (seek($self->{'_FH'}, 0, 0)) {
		carp "Error on reading file";
		return 0;
	}
	# Read the first four bytes of the file
	# and interpret them as ASCII characters
	read($self->{'_FH'}, $file_signature, 4) or return 0;
	$file_signature = unpack('A4', $file_signature);
	if ($file_signature eq 'FIBA') {
		print STDERR "Probably, an ABIF file stored in little endian order\n";
		print STDERR "Unsupported ABIF file structure (because deprecated)\n";
	}
	return ($file_signature eq 'ABIF');
}

=head2 abif_version()

  Usage    : my $v = $ab1->abif_version();
  Returns  : The ABIF file version number (e.g., '1.01').

Used to determine the ABIF file version number.
  
=cut

sub abif_version {
	my $self = shift;
	my $version;
	unless (defined $self->{'_ABIF_VERSION'}) {
		# Version number is stored in bytes 4 and 5
		seek($self->{'_FH'}, 4, 0) or croak "Error on reading file";
		read($self->{'_FH'}, $version, 2) or croak "Error on reading file";
		$version = unpack('n', $version);
		$self->{'_ABIF_VERSION'} = $version / 100;
	}
	return $self->{'_ABIF_VERSION'};
}

=head2 num_dir_entries()

  Usage    : my $n = $ab1->num_dir_entries();
  Returns  : The number of data items contained in the file.
  
Used to determine the number of directory entries in the ABIF file.

=cut

sub num_dir_entries {
	my $self = shift;
	return $self->{'_NUMELEM'};
}

=head2 data_offset()

  Usage    : my $n = $ab1->data_offset();
  Returns  : The offset, in bytes, of the first directory entry
             with respect to the beginning of the file
  
Used to determine the data offset of the directory entries.

=cut

sub data_offset {
	my $self = shift;
	return $self->{'_DATAOFFSET'};
}

=head2 get_directory()

  Usage    : my %DirEntry = $ab1->get_directory($tag_name, $tag_num);
  Returns  : a hash with the content of the specified directory entry,
             or () if the tag is not found.
             
Retrieves the directory entry identified by the pair ($tag_name, $tag_num).
The $tag_name is a four letter ASCII code and $tag_num is an integer
(typically, 1 <= $tag_num <= 1000). The returned hash has the following keys:

                TAG_NAME: the tag name;
                TAG_NUMBER: the tag number;
                ELEMENT_TYPE: a string indicating the type of the data item
                              ('char', 'byte', 'float', 'pString', etc...);
                ELEMENT_SIZE: the size, in bytes, of one element;
                NUM_ELEMENTS: the number of elements in the data item;
                DATA_SIZE: the size, in bytes, of the data item;
                DATA_ITEM: the raw sequence of bytes of the data item.

Nota Bene: it is upon the caller to interpret the data item field correctly
(typically, by unpack()ing the item).

Refer to the L</"SEE ALSO"> Section for further information.

=cut

sub get_directory {
	my ($self, $tag_name, $tag_number) = @_;
	my %DirEntry;
	my $field;
	my $raw_data;

	if ($self->search_tag($tag_name, $tag_number)) { # Found!
		$DirEntry{TAG_NAME} = $tag_name;
		$DirEntry{TAG_NUMBER} = $tag_number;
		# Read and unpack the remaining bytes
		read($self->{'_FH'}, $field, 12);
		#my ($et, $es, $ne, $ds) = unpack('ssll', $field);
		my ($et, $es, $ne, $ds) = unpack('nnNN', $field);
		# Element type code (signed 16 bit integer)		
		if ($et > 1023) {
			$DirEntry{ELEMENT_TYPE} = 'user';
		}
		else {
			$DirEntry{ELEMENT_TYPE} = $self->{TYPES}{$et};
		}
		# Element size (signed 16 bit integer)
		$DirEntry{ELEMENT_SIZE} = $es;
		# Number of element in this item (signed 32 bit integer)
		$DirEntry{NUM_ELEMENTS} = $ne;
		# Size of the item in bytes (signed 32 bit integer)
		$DirEntry{DATA_SIZE} = $ds;
		# Get data item
		if ($DirEntry{DATA_SIZE} > 4) {
		 	# The data item position is given by the data offset field
		 	read($self->{'_FH'}, $field, 4);
		 	#$field = unpack('l', $field);
		 	$field = unpack('N', $field);
			seek($self->{'_FH'}, $field, 0);
			read($self->{'_FH'}, $raw_data, $DirEntry{DATA_SIZE});
		}
		else {
			# if data size <= 4 then the data item is stored in the data offset field itself
			# (the current file handle position)
			read($self->{'_FH'}, $raw_data, $DirEntry{DATA_SIZE});
		}
		$DirEntry{DATA_ITEM} = $raw_data; # Return raw data

		return %DirEntry;
	}
	
	return ();
}

=head2 get_data_item()

  Usage    : my @data = $ab1->get_data_item($tag_name, $tag_num, $template);
  Returns  : (), if the tag is not found; otherwise:
             a list of elements unpacked according to $template.
  
Retrieves the data item specified by the pair ($tag_name, $tag_num) and,
unpacks it according to $template. The $tag_name is a four letter
ASCII code and $tag_num is an integer (typically, 1 <= $tag_num <= 1000).
The $template has the same format as in the pack() function.               
  
Refer to the L</"SEE ALSO"> Section for further information.

=cut

sub get_data_item {
	my $self = shift;
	my $tag_name = shift;
	my $tag_number = shift;
	my $template = shift;
	my $field;
	my $raw_data;

	if ($self->search_tag($tag_name, $tag_number)) { # Found!
		# Read the remaining bytes of the current directory entry
		read($self->{'_FH'}, $field, 12);
		# Unpack data size
		#my ($data_size) = unpack('x8l', $field);
		my ($data_size) = unpack('x8N', $field);
		if ($data_size > 4) {
		 	# The data item position is given by the data offset field
		 	read($self->{'_FH'}, $field, 4);
		 	#$field = unpack('l', $field);
		 	$field = unpack('N', $field);
			seek($self->{'_FH'}, $field, 0);
			read($self->{'_FH'}, $raw_data, $data_size);
		}
		else {
			# if data size <= 4 then the data item is stored in the data offset field itself
			# (the current file handle position)
			read($self->{'_FH'}, $raw_data, $data_size);			
		}
		my @data = unpack($template, $raw_data);
		return @data;
	}	
	return ();
}

=head2 search_tag()

  Usage    : if (search_tag($tag_name, $tag_num)) { # etc...
  Returns  : 1 if the tag is found; 0, otherwise

Performs a linear scan of the directory entries until the specified data tag
is matched.

If the tag is found then the file handle is positioned 
just after thet tag number (ready to read the element type).

=cut

sub search_tag {
	my ($self, $tag_name, $tag_number) = @_;
	my ($t1, $t2, $field);
	my $offset = ${$self->{'_TAG_INDEX'}}{$tag_name . $tag_number};
	if (defined $offset) {
		seek($self->{'_FH'}, $offset + 8, 0);
		return 1;
	}
	else { # Linear search
		my $i = 0;
		do {
			seek($self->{'_FH'}, $self->{'_DATAOFFSET'} + ($DIR_SIZE * $i), 0);
			# Read Tag name and number (8 bytes)
			read($self->{'_FH'}, $field, 8);
			#($t1, $t2) = unpack('A4l', $field);
			($t1, $t2) = unpack('A4N', $field); # Four ASCII letters and a signed 32 bit integer
			$i++;
		}
		while ((($t1 ne $tag_name) or ($t2 != $tag_number)) and ($i < $self->{'_NUMELEM'}));  

		return ($i < $self->{'_NUMELEM'});
	}
}

=head1 AB1 COMMON TAGS

The following methods work for files from the AB 3730/3730xl
Data Collection Software v2.0 and v3.0 on the Applied Biosystems 3730/3730xl
Genetic/DNA Analyzer or from ABI Prism(R) 3100/3100-Avant Analyzer
Data Collection Software v2.0.

=head2 sample_name()

  Usage    : my $name = sample_name();
  Returns  : a string containing the sample name;
             '' if this tag is not in the file.
  
=cut

sub sample_name {
	my $self = shift;
	unless (defined $self->{'_SMPL1'}) {
		my ($sn) = $self->get_data_item('SMPL', 1, 'xA*');
		$self->{'_SMPL1'} = ($sn) ? $sn : ''; 
	}
	return $self->{'_SMPL1'};
}

=head2 data_collection_software_version()

  Usage    : my $v = data_collection_software_version();
  Returns  : the data collection software version.
             '' if this tag is not in the file.

=cut

sub data_collection_software_version {
	my $self = shift;
	unless (defined $self->{'_SVER1'}) {
		my ($sw) = $self->get_data_item('SVER', 1, 'xA*');
		$self->{'_SVER1'} = ($sw) ? $sw : '';
	}
	return $self->{'_SVER1'};
}

=head2 data_collection_firmware_version()

  Usage    : my $v = data_collection_firmware_version();
  Returns  : the data collection firmware version;
             '' if this tag is not in the file.

=cut

sub data_collection_firmware_version {
	my $self = shift;
	unless (defined $self->{'_SVER3'}) {
		my ($fw) = $self->get_data_item('SVER', 3, 'xA*');
		$self->{'_SVER3'} = ($fw) ? $fw : '';
	}
	return $self->{'_SVER3'};
}

=head2 official_instrument_name()

  Usage    : my $v = official_instrument_name();
  Returns  : the official instrument name;
             '' if this tag is not in the file.


=cut

sub official_instrument_name {
	my $self = shift;
	unless (defined $self->{'_HCFG3'}) {
		my ($name) = $self->get_data_item('HCFG', 3, 'xA*');
		$self->{'_HCFG3'} = ($name) ? $name : '';
	}
	return $self->{'_HCFG3'};	
}

=head2 well_id()

  Usage    : my $well_id = well_id();
  Returns  : the well ID;
             '' if this tag is not in the file.
  
=cut

sub well_id {
	my $self = shift;
	unless (defined $self->{'_TUBE1'}) {
		my ($wid) = $self->get_data_item('TUBE', 1, 'xA*');
		$self->{'_TUBE1'} = ($wid) ? $wid : '';
	}
	return $self->{'_TUBE1'};
}

=head2 capillary_number()

  Usage    : my $cap_n = capillary_number();
  Returns  : the LANE/Capillary number;
             0 if this tag is not in the file.
  
=cut

sub capillary_number {
	my $self = shift;
	unless (defined $self->{'_LANE1'}) {
		my ($cap_n) = $self->get_data_item('LANE', 1, 'n');
		$self->{'_LANE1'} = ($cap_n) ? $cap_n : 0;
	}
	return $self->{'_LANE1'};
}


=head2 user()

  Usage    : my $user = user();
  Returns  : the name of the user who created the plate;
             '' if this tag is not in the file.
  
Nota Bene: this field is optional

=cut

sub user {
	my $self = shift;
	unless (defined $self->{'_User1'}) {
		my ($u) = $self->get_data_item('User', 1, 'xA*');
		$self->{'_User1'} = ($u) ? $u : '';
	}
	return $self->{'_User1'};
}

=head2 instrument_name_and_serial_number()

  Usage    : my $sn = instrument_name_and_serial_number()
  Returns  : a string with the instrument name and the serial number;
             '' if this tag is not in the file.
 
=cut

sub instrument_name_and_serial_number {
	my $self = shift;
	unless (defined $self->{'_MCHN1'}) {
		my ($sn) = $self->get_data_item('MCHN', 1, 'xA*');
		$self->{'_MCHN1'} = ($sn) ? $sn : '';
	}
	return $self->{'_MCHN1'};
}

=head2 sequencing_analysis_param_filename()

  Usage    : my $f = sequencing_analysis_param_filename();
  Returns  : the Sequencing Analysis parameters filename;
             '' if this tag is not in the file.
  
=cut

sub sequencing_analysis_param_filename {
	my $self = shift;
	unless (defined $self->{'_APFN2'}) {
		my ($sa) = $self->get_data_item('APFN', 2, 'xA*');
		$self->{'_APFN2'} = ($sa) ? $sa : '';
	}
	return $self->{'_APFN2'};
}

=head2 comment()

  Usage    : my $comment = comment();
  Returns  : the comment associated to the file;
             '' if this tag is not in the file.
  
This is an optional field.

=cut

sub comment {
	my $self = shift;
	unless (defined $self->{'_CMNT1'}) {
		my ($c) = $self->get_data_item('CMNT', 1, 'xA*');
		$self->{'_CMNT1'} = ($c) ? $c : '';
	}
	return $self->{'_CMNT1'};
}

=head2 base_order()

  Usage    : my @bo = $ab->base_order();
  Returns  : the order in which the bases are stored in the file;
             () if this tag is not in the file.
  
=cut

sub base_order {
	my $self = shift;
	unless (defined $self->{'_FWO_1'}) {
		my ($bases) = $self->get_data_item('FWO_', 1, 'A*');
		if (defined $bases) {
			my @bo = split('', $bases);
			$self->{'_FWO_1'} = [ @bo ];			
		}
		else {
			$self->{'_FWO_1'} = [ ];
		}
	}
	return @{$self->{'_FWO_1'}};
}

=head2 order_base()

  Usage    : my %bases = $ab->order_base();
  Returns  : the indices of the bases, as they stored in the file;
             () if the base order is not present in the file.
  
This method does the opposite as base_order() does. 
  
=cut

sub order_base {
	my $self = shift;
	unless (defined $self->{'_OB'}) {
		my @bo = $self->base_order();
		my %ob = ();
		for (my $i = 0; $i < scalar(@bo); $i++) {
			$ob{$bo[$i]} = $i;
		}
		$self->{'_OB'} = ( %ob );
	}
	return $self->{'_OB'};
}

=head2 num_capillaries()

  Usage    : my $nc = num_capillaries();
  Returns  : the number of capillaries;
             0 if this tag is not in the file.
  
=cut

sub num_capillaries {
	my $self = shift;
	unless (defined $self->{'_NLNE1'}) {
		my ($nc) = $self->get_data_item('NLNE', 1, 'n');
		$self->{'_NLNE1'} = ($nc) ? $nc : 0;
	}
	return $self->{'_NLNE1'};
}

=head2 sample_tracking_id()

  Usage    : my $sample_id = sample_tracking_id();
  Returns  : the sample tracking ID;
             '' if this tag is not in the file.
  
=cut

sub sample_tracking_id {
	my $self = shift;
	unless (defined $self->{'_LIMS1'}) {
		my ($sid) = $self->get_data_item('LIMS', 1, 'xA*');
		$self->{'_LIMS1'} = ($sid) ? $sid : '';
	}
	return $self->{'_LIMS1'};
}

=head2 analysis_protocol_xml()

  Usage    : my $xml = analysis_protocol_xml();
  Returns  : the Analysis Protocol XML string;
             '' if this tag is not in the file.
  
=cut

sub analysis_protocol_xml {
	my $self = shift;
	unless (defined $self->{'_APrX1'}) {
		my ($apxml) = $self->get_data_item('APrX', 1, 'c*');
		$self->{'_APrX1'} = ($apxml) ? $apxml : '';
	}
	return $self->{'_APrX1'};
}


=head2 analysis_protocol_settings_name()

  Usage    : my $name = analysis_protocol_settings_name();
  Returns  : the Analysis Protocol settings name;
             '' if this tag is not in the file.
  
=cut

sub analysis_protocol_settings_name {
	my $self = shift;
	unless (defined $self->{'_APrN1'}) {
		my ($apsn) = $self->get_data_item('APrN', 1, 'Z*');
		$self->{'_APrN1'} = ($apsn) ? $apsn : '';
	}
	return $self->{'_APrN1'};
}

=head2 analysis_protocol_settings_version()

  Usage    : my $name = analysis_protocol_settings_version();
  Returns  : the Analysis Protocol settings version;
             '' if this tag is not in the file.
  
=cut

sub analysis_protocol_settings_version {
	my $self = shift;
	unless (defined $self->{'_APrV1'}) {
		my ($apsv) = $self->get_data_item('APrV', 1, 'Z*');
		$self->{'_APrV1'} = ($apsv) ? $apsv : '';
	}
	return $self->{'_APrV1'};
}

=head2 analysis_protocol_xml_schema_version()

  Usage    : my $name = analysis_protocol_xml_schema_version();
  Returns  : the Analysis Protocol XML schema version;
             '' if this tag is not in the file.
  
=cut

sub analysis_protocol_xml_schema_version {
	my $self = shift;
	unless (defined $self->{'_APXV1'}) {
		my ($apxv) = $self->get_data_item('APXV', 1, 'Z*');
		$self->{'_APXV1'} = ($apxv) ? $apxv : '';
	}
	return $self->{'_APXV1'};
}

=head2 results_group()

  Usage    : my $name = results_group();
  Returns  : the results group name;
             '' if this tag is not in the file.
  
=cut

sub results_group {
	my $self = shift;
	unless (defined $self->{'_RGNm1'}) {
		my ($rg) = $self->get_data_item('RGNm', 1, 'Z*');
		$self->{'_RGNm1'} = ($rg) ? $rg : '';
	}
	return $self->{'_RGNm1'};
}

=head2 run_module_name()

  Usage    : my $name = run_module_name();
  Returns  : the run module name;
             '' if this tag is not in the file.
  
=cut

sub run_module_name {
	my $self = shift;
	unless (defined $self->{'_RMdN1'}) {
		my ($rn) = $self->get_data_item('RMdN', 1, 'Z*');
		$self->{'_RMdN1'} = ($rn) ? $rn : '';
	}
	return $self->{'_RMdN1'};	
}

=head2 run_module_version()

  Usage    : my $name = run_module_version();
  Returns  : the run module version;
             '' if this tag is not in the file.
  
=cut

sub run_module_version {
	my $self = shift;
	unless (defined $self->{'_RMdV1'}) {
		my ($rv) = $self->get_data_item('RMdV', 1, 'Z*');
		$self->{'_RMdV1'} = ($rv) ? $rv : '';
	}
	return $self->{'_RMdV1'};	
}

=head2 raw_data_for_channel()

  Usage    : my @data = raw_data_for_channel($channel_number);
  Returns  : the channel $channel_number raw data;
             () if the channel number is out of range
             or the tag is not in the file.
  
There are four channels in an ABIF file, numbered from 1 to 4.
An optional channel number 5 exists in some files. If a channel
is specified out of such ranges, the empty list is returned.
If channel 5 data is requested, but such data is not found, then
returns undef.

=cut

sub  raw_data_for_channel {
	my ($self, $channel_number) = @_;
	my $k = '_DATA' . $channel_number;
	unless (defined $self->{$k}) {
		if ($channel_number > 0 and $channel_number <= 5) {
			$channel_number = 105 if ($channel_number == 5);
			my @data = $self->get_data_item('DATA', $channel_number, 'n*');
			$self->{$k} = (@data) ? [ @data ] : [ ];
		}
		else {
			$self->{$k} = [ ];
		}
	}
	return @{$self->{$k}};
}

=head2 num_dyes()

  Usage    : my $n = num_dyes();
  Returns  : the number of dyes;
             0 if this tag is not in the file.
  
=cut
sub num_dyes {
	my $self = shift;
	unless (defined $self->{'_Dye#1'}) {
		my ($nd) = $self->get_data_item('Dye#', 1, 'n');
		$self->{'_Dye#1'} = ($nd) ? $nd : 0;
	}
	return $self->{'_Dye#1'};
}

=head2 dye_name()

  Usage    : my $n = dye_name($n);
  Returns  : the name of dye number $n;
             '' if this tag is not in the file or $n is not
             in the range [1..4].
             
=cut

sub dye_name {
	my ($self, $n) = @_;
	my $k = '_DyeN'. $n;
	unless (defined $self->{$k}) {
		if ($n > 0 and $n < 5) {
			my ($dn) = $self->get_data_item('DyeN', $n, 'xA*');
			$self->{$k} = ($dn) ? $dn : '';
		}
		else {
			$self->{$k} = '';
		}
	}
	return $self->{$k};
}

=head1 AB1 NEWER TAGS

The following methods work for files from the AB 3730/3730xl
Data Collection Software v3.0 on the Applied Biosystems 3730/3730xl
Genetic Analyzer.

=head2 container_owner()

  Usage    : my $owner = container_owner();
  Returns  : the container's owner;
             '' if this tag is not in the file.
  
=cut

sub container_owner {
	my $self = shift;
	unless (defined $self->{'_CTOw1'}) {
		my ($co) = $self->get_data_item('CTOw', 1, 'Z*');
		$self->{'_CTOw1'} = ($co) ? $co : '';
	}
	return $self->{'_CTOw1'};
}

=head1 SEQSCAPE V2.5 AND SEQUENCING ANALYSIS V5.2 TAGS

The following methods work from file processed by SeqScape v2.5
or Sequencing Analysis v5.2.

=head2 quality_values()

  Usage    : my @qv = quality_values();
  Returns  : the list of quality values;
             () if this tag is not in the file.
  
=cut

sub quality_values {
	my $self = shift;
	unless (defined $self->{'_PCON2'}) {
		# Load and cache quality values
		my @qv = $self->get_data_item('PCON', 2, 'c*');
		$self->{'_PCON2'} = (@qv) ? [ @qv ] : [ ];
	}
	return @{$self->{'_PCON2'}};
}

=head2 quality_values_ref()

  Usage    : my $ref_to_qv = quality_values_ref();
  Returns  : a reference to the list of quality values;
             a reference to the empty list if this tag is not in the file.
  
=cut

sub quality_values_ref {
	my $self = shift;
	unless (defined $self->{'_PCON2'}) {
		# Load and cache quality values
		my @qv = $self->get_data_item('PCON', 2, 'c*');
		$self->{'_PCON2'} = (@qv) ? [ @qv ] : [ ];
	}
	return $self->{'_PCON2'};
}

=head2 edited_quality_values()

  Usage    : my @qv = edited_quality_values();
  Returns  : the list of edited quality values;
             () if this tag is not in the file.
  
=cut

sub edited_quality_values {
	my $self = shift;
	unless (defined $self->{'_PCON1'}) {
		my @qv = $self->get_data_item('PCON', 1, 'c*');
		$self->{'_PCON1'} = (@qv) ? [ @qv ] : [ ];
	}
	return @{$self->{'_PCON2'}};
}

=head2 edited_quality_values_ref()

  Usage    : my $ref_to_qv = edited_quality_values_ref();
  Returns  : a reference to the list of edited quality values;
             a reference to the empty list if this tag is not in the file.
  
=cut

sub edited_quality_values_ref {
	my $self = shift;
	unless (defined $self->{'_PCON1'}) {
		my @qv = $self->get_data_item('PCON', 1, 'c*');
		$self->{'_PCON1'} = (@qv) ? [ @qv ] : [ ];
	}
	return $self->{'_PCON2'};
}


=head2 sequence()

  Usage    : my $sequence = sequence();
  Returns  : the string of the base called sequence;
             '' if this tag is not in the file.
  
=cut

sub sequence {
	my $self = shift;
	unless (defined $self->{'_PBAS2'}) {
		my ($seq) = $self->get_data_item('PBAS', 2, 'A*');
		$self->{'_PBAS2'} = ($seq) ? $seq : '';
	}
	return $self->{'_PBAS2'};
}

=head2 sequence_length()

  Usage    : my $l = sequence_length();
  Returns  : the length of the base called sequence;
             0 if the sequence is not in the file.
  
=cut

sub sequence_length() {
	my $self = shift;
	my $seq = $self->sequence();
	return 0 unless defined $seq;
	return length($seq);
}

=head2 edited_sequence()

  Usage    : my $sequence = edited_sequence();
  Returns  : the string of the edited base called sequence;
             '' if this tag is not in the file.
  
=cut

sub edited_sequence {
	my $self = shift;
	unless (defined $self->{'_PBAS1'}) {
		my ($seq) = $self->get_data_item('PBAS', 1, 'A*');
		$self->{'_PBAS1'} = ($seq) ? $seq : '';
	}
	return $self->{'_PBAS1'};
}

=head2 edited_sequence_length()

  Usage    : my $l = edited_sequence_length();
  Returns  : the length of the base called sequence;
             0 if the sequence is not in the file.
  
=cut

sub edited_sequence_length() {
	my $self = shift;
	my $seq = $self->edited_sequence();
	return 0 unless defined $seq;
	return length($seq);
}


=head2 analyzed_data_for_channel()

  Usage    : my @data = analyzed_data_for_channel($channel_number);
  Returns  : the channel $channel_number analyzed data;
             () if the channel number is out of range
             or this tag is not in the file.
  
There are four channels in an ABIF file, numbered from 1 to 4.
An optional channel number 5 exists in some files.
If a channel is specified out of such ranges, the empty list is returned.
If channel 5 data is requested, but such data is not found, returns undef.

=cut

sub  analyzed_data_for_channel {
	my ($self, $channel_number) = @_;
	my $key = '_DATA' . $channel_number;
	unless (defined $self->{$key}) {
		if ($channel_number > 0 and $channel_number <= 5) {
			$channel_number = 197 if ($channel_number == 5);
			my @data = $self->get_data_item('DATA', $channel_number + 8, 'n*');
			$self->{$key} = (@data) ? [ @data ] : [ ];
		}
		else {
			$self->{$key} = [ ];
		}
	}
	return @{$self->{$key}};
}


=head2 peak1_location_orig()

  Usage    : my $pl = peak1_location_orig();
  Returns  : The peak 1 location (orig);
             0 if this tag is not in the file.
  
=cut

sub peak1_location_orig {
	my $self = shift;
	unless (defined $self->{'_B1Pt1'}) {
		my ($pl) = $self->get_data_item('B1Pt', 1, 'n');
		$self->{'_B1Pt1'} = ($pl) ? $pl : 0;
	}
	return $self->{'_B1Pt1'};
}

=head2 peak1_location()

  Usage    : my $pl = peak1_location();
  Returns  : The peak 1 location;
             0 if this tag is not in the file.
  
=cut

sub peak1_location {
	my $self = shift;
	unless (defined $self->{'_B1Pt2'}) {
		my ($pl) = $self->get_data_item('B1Pt', 2, 'n');
		$self->{'_B1Pt2'} = ($pl) ? $pl : 0;
	}
	return $self->{'_B1Pt2'};
}

=head2 base_spacing()

  Usage    : my $spacing = base_spacing();
  Returns  : the spacing (a float);
             0.0 if this tag is not in the file.

=cut

sub base_spacing() {
	my $self = shift;
	unless (defined $self->{'_SPAC3'}) {
		my ($s) = $self->get_data_item('SPAC', 3, 'B32');
		if ($s) {
			$self->{'_SPAC3'} = $self->_ieee_single_prec_float($s);
		}
		else {
			$self->{'_SPAC3'} = 0.0;
		}
	}
	return $self->{'_SPAC3'};
}

=head2 basecaller_version()

  Usage    : my $v = basecaller_version();
  Returns  : a string indicating the basecaller version (e.g., 'KB 1.3.0');
             '' if this tag is not in the file.
 
=cut

sub basecaller_version {
	my $self = shift;
	unless (defined $self->{'_SVER2'}) {
		# Basecaller version is a pString (Pascal string):
		# the first byte (which is skipped) is the length of the string
		my ($v) = $self->get_data_item('SVER', 2, 'xA*');
		$self->{'_SVER2'} = ($v) ? $v : '';
	}
	return $self->{'_SVER2'};
}

=head2 base_locations()

  Usage    : my @bl = base_locations();
  Returns  : the list of base locations;
             () if this tag is not in the file.
  
=cut

sub base_locations {
	my $self = shift;
	unless (defined $self->{'_PLOC2'}) {
		my @bl = $self->get_data_item('PLOC', 2, 'n*');
		$self->{'_PLOC2'} = (@bl) ? [ @bl ] : [ ];
	}
	return @{$self->{'_PLOC2'}};
}

=head2 base_locations_edited()

  Usage    : my @bl = base_locations_edited();
  Returns  : the list of base locations (edited);
             () if this tag is not in the file.

=cut

sub base_locations_edited {
	my $self = shift;
	unless (defined $self->{'_PLOC1'}) {
		my @bl = $self->get_data_item('PLOC', 1, 'n*');
		$self->{'_PLOC1'} = (@bl) ? [ @bl ] : [ ];
	}
	return @{$self->{'_PLOC1'}};
}

=head2 basecaller_bcp_dll()

  Usage    : my $v = basecaller_bcp_dll();
  Returns  : a string with the basecalled BCP/DLL;
             '' if this tag is not in the file.
  
=cut

sub basecaller_bcp_dll {
	my $self = shift;
	unless (defined $self->{'_SPAC2'}) {
		my ($v) = $self->get_data_item('SPAC', 2, 'xA*');
		$self->{'_SPAC2'} = ($v) ? $v : '';
	}
	return $self->{'_SPAC2'};	
}

=head2 signal_level()

  Usage    : my %signal_level = signal_level();
  Returns  : the signal level for each dye;
             () if this tag is not in the file.
             
=cut

sub signal_level {
	my $self = shift;
	unless (defined $self->{'_S/N%1'}) {
		my @sl = $self->get_data_item('S/N%', 1, 'n*');
		unless (@sl) {
			$self->{'_S/N%1'} = { };
		}
		else {
			my %signal = ();
			my @bo = $self->base_order();
			for (my $i = 0; $i < scalar(@sl); $i++) {
				$signal{$bo[$i]} = $sl[$i];
			}
			$self->{'_S/N%1'} = { %signal };
		}
	}
	return %{$self->{'_S/N%1'}};
}

=head2 raw_trace()

  Usage    : my @trace = raw_trace($base);
  Returns  : the raw trace corresponding to base $base;
             () if the data is not in the file.
  
The possible values for $base are 'A', 'C', 'G' and 'T' (case insensitive).

=cut

sub raw_trace {
	my ($self, $base) = @_;
	
	$base =~ /^[ACGTacgt]$/ or return ();
	$base = uc($base);
	return $self->raw_data_for_channel($self->order_base($base) + 1);
}

=head2 trace()

  Usage    : my @trace = trace($base);
  Returns  : the (analyzed) trace corresponding to base $base;
             () if the data is not in the file.
  
The possible values for $base are 'A', 'C', 'G' and 'T'.

=cut

sub trace {
	my ($self, $base) = @_;
	
	$base =~ /^[ACGTacgt]$/ or return ();
	$base = uc($base);
	return $self->analyzed_data_for_channel($self->order_base($base) + 1);
}

=head2 noise()

  Usage    : my %noise = $ab->noise();
  Returns  : the estimated noise for each dye;
             () if this tag is not in the file.
  
This method works only with files containing data processed by the KB Basecaller.

=cut

sub noise {
	my $self = shift;
	unless (defined $self->{'_NOIS1'}) {
		my %noise = ();	
		my ($bits) = $self->get_data_item('NOIS', 1, 'B*');
		unless (defined $bits) {
			$self->{'_NOIS1'} = { };
		}
		else {
			my @bo = $self->base_order();
			for (my $i = 0; $i < length($bits); $i += 32) {
				# Convert to float
				$noise{$bo[$i / 32]} = $self->_ieee_single_prec_float(substr($bits, $i, 32));
			}
			$self->{'_NOIS1'} = { %noise };
		}
	}	
	return %{$self->{'_NOIS1'}};
}

=head1 METHODS FOR ASSESSING QUALITY

The following methods compute some values that help assessing the quality
of the data.

=head2 avg_signal_to_noise_ratio()

  Usage    : my $sn_ratio = $ab->avg_signal_to_noise_ratio()
  Returns  : the average signal to noise ratio (only for the KB Basecaller);
             0 if the information needed to compute such value is missing.

This method works only with files containing data processed by the KB Basecaller.

=cut

sub avg_signal_to_noise_ratio {
	my $self = shift;
	my %sl = $self->signal_level();
	return 0 unless %sl;
	my %noise = $self->noise();
	return 0 unless %noise;
	my $avg = 0;
	foreach my $base (keys %sl) {
		$avg += $sl{$base} / $noise{$base};
	}
	return $avg / scalar(keys %sl);
}

=head2 clear_range_start()

  Usage    : my $cl_start = $ab->clear_range_start();
             my $cl_start = $ab->clear_range_start($window_width,
                                                   $bad_bases_threshold,
                                                   $quality_threshold);
  Returns  : the clear range start position (counting from zero);
             -1 if $window_width is greater than the number of quality values;
             -1 if the information needed to compute such value is missing.
  
The Sequencing Analysis program determines the clear range of the sequence 
by trimming bases from the 5' to 3' ends until fewer than 4 bases out of 20
have a quality value less than 20. You can change these parameters by explicitly
passing arguments to this method (the default values are $window_width = 20,
$bad_bases_threshold = 4, $quality_threshold = 20). Note that Sequencing Analysis
counts the bases starting from one, so you have to add one to the return values
to get consistent results.

=cut

sub clear_range_start {
	my $self = shift;
	my $window = 20;
	my $bad_bases = 4;
	my $threshold = 20;
	if (@_) {
		$window = shift;
		$bad_bases = shift;
		$threshold = shift;
	}
	my $qv_ref = $self->quality_values_ref();
	return -1 unless defined $qv_ref;

	my $N = scalar(@$qv_ref);
	return -1 if ($N < $window); # Not enough quality values
	
	my $j;
	my $n = 0; # of bases with qv < 20
	for ($j = 0; $j < $window; $j++) {
		if ($$qv_ref[$j] < $threshold) {
			$n++;
		}
	}
	return 0 if ($n < $bad_bases);
	
	do {
		if ($$qv_ref[$j - $window] < $threshold) {
			$n--;
		}
		if ($$qv_ref[$j] < $threshold) {
			$n++;
		}
		$j++;
	}
	while ($n >= $bad_bases and $j < $N);
	return $j - $window;
}


=head2 clear_range_stop()

  Usage    : my $cl_stop = $ab->clear_range_stop();
             my $cl_stop = $ab->clear_range_stop($window_width,
                                                 $bad_bases_threshold,
                                                 $quality_threshold);
  Returns  : the clear range stop position (counting from zero);
             -1 if $window_width is greater than the number of quality values;
             -1 if the information needed to compute such value is missing.
  
The Sequencing Analysis program determines the clear range of the sequence 
by trimming bases from the 5' to 3' ends until fewer than 4 bases out of 20
have a quality value less than 20. You can change these parameters by explicitly
passing arguments to this method (the default values are $window_width = 20,
$bad_bases_threshold = 4, $quality_threshold = 20).

=cut

sub clear_range_stop {
	my $self = shift;
	my $window = 20;
	my $bad_bases = 4;
	my $threshold = 20;
	if (@_) {
		$window = shift;
		$bad_bases = shift;
		$threshold = shift;
	}
	my $qv_ref = $self->quality_values_ref();
	return -1 unless  defined $qv_ref;
	
	my $N = scalar(@$qv_ref);
	return -1 if ($N < $window); # Not enough quality values
	
	my $j;
	my $n = 0; # of bases with qv < 20
	for ($j = $N - 1; $j >= $N - $window; $j--) {
		if ($$qv_ref[$j] < $threshold) {
			$n++;
		}
	}
	return $j + $window if ($n < $bad_bases);
	
	do {
		if ($$qv_ref[$j + $window] < $threshold) {
			$n--;
		}
		if ($$qv_ref[$j] < $threshold) {
			$n++;
		}
		$j--;
	}
	while ($n >= $bad_bases and $j >= 0);
	return $j + $window;
}


=head2 sample_score()

  Usage    : my $ss = $ab->sample_score();
           : my $ss = $ab->sample_score($window_width, $bad_bases_threshold,
                                        $quality_threshold);
  Returns  : the sample score associated to the sequence;
             0 if the information needed to compute such value is missing.
  
The sample score is the average quality value of the bases
in the clear range of the sequence. See the clear_range_start() and
clear_range_stop() methods for the meaning of the optional arguments of
this method.

=cut

sub sample_score {
	my $self = shift;
	my $start;
	my $stop;
	if (@_) {
		my $window = shift;
		my $bad_bases = shift;
		my $threshold = shift;
		$start = $self->clear_range_start($window, $bad_bases, $threshold);
		$stop  = $self->clear_range_stop($window, $bad_bases, $threshold);
	}
	else {
		$start = $self->clear_range_start();
		$stop = $self->clear_range_stop();
	}
	my $qv_ref = $self->quality_values_ref();
	return 0 unless ($start >= 0) and ($stop >= 0) and defined $qv_ref;
	# Compute average quality value in the clear range
	my $sum = 0;
	for (my $i = $start; $i <= $stop; $i++) {
		$sum += $$qv_ref[$i];
	}
	return $sum / ($stop - $start + 1);
}

=head2 num_medium_quality_bases()

  Usage    : my $n = $ab->num_medium_quality_bases($min_qv, $max_qv, $start, $stop);
             my $n = $ab->num_medium_quality_bases($min_qv, $max_qv);
  Returns  : the number of bases in the range [$start, $stop], or in the whole
             sequence if no range is specified, with $min_qv <= quality value <= $max_qc;
             -1 if the information needed to compute such value is missing.

=cut

sub num_medium_quality_bases {
	my $self = shift;
	my $min = shift;
	my $max = shift;
	my $start;
	my $stop;
	if (@_) {
		$start = shift;
		$stop = shift;
	}
	else {
		$start = 0;
		$stop = -1;
	}

	my $qv_ref = $self->quality_values_ref();
	return -1 unless defined $qv_ref;
	my $n = 0;
	if ($start <= $stop) {
		for (my $i = $start; $i <= $stop; $i++) {
			if ($$qv_ref[$i] >= $min and $$qv_ref[$i] <= $max) {
				$n++;
			}
		}
	}
	else { # Count all the quality values
		foreach my $qv (@$qv_ref) {
			if ($qv >= $min and $qv <= $max) {
				$n++;
			}
		}
	}
	return $n;
}

=head2 num_high_quality_bases()

  Usage    : my $n = $ab->num_high_quality_bases($threshold, $start, $stop);
             my $n = $ab->num_high_quality_bases($threshold);
  Returns  : the number of bases in the range [$start, $stop], or in the whole
             sequence if no range is specified, with quality value >= $threshold;
             -1 if the information needed to compute such value is missing.
             
=cut

sub num_high_quality_bases {
	my $self = shift;
	my $min = shift;
	my $start;
	my $stop;
	if (@_) {
		$start = shift;
		$stop = shift;
	}
	else {
		$start = 0;
		$stop = -1;
	}
	
	my $qv_ref = $self->quality_values_ref();
	return -1 unless defined $qv_ref;
	my $n = 0;
	if ($start <= $stop) {
		for (my $i = $start; $i <= $stop; $i++) {
			if ($$qv_ref[$i] >= $min) {
				$n++;
			}
		}
	}
	else { # Count all the quality values
		foreach my $qv (@$qv_ref) {
			if ($qv >= $min) {
				$n++;
			}
		}
	}
	return $n;
}

=head2 num_low_quality_bases()

  Usage    : my $n = $ab->num_low_quality_bases($threshold, $start, $stop);
             my $n = $ab->num_low_quality_bases($threshold);
  Returns  : the number of bases in the range [$start, $stop], or in the whole
             sequence if no range is specified, with quality value <= $threshold;
             -1 if the information needed to compute such value is missing.
             
=cut

sub num_low_quality_bases {
	my $self = shift;
	my $max = shift;
	my $start;
	my $stop;
	if (@_) {
		$start = shift;
		$stop = shift;
	}
	else {
		$start = 0;
		$stop = -1;
	}
	
	my $qv_ref = $self->quality_values_ref();
	return -1 unless defined $qv_ref;
	my $n = 0;
	if ($start <= $stop) {
		for (my $i = $start; $i <= $stop; $i++) {
			if ($$qv_ref[$i] <= $max) {
				$n++;
			}
		}
	}
	else { # Count all the quality values
		foreach my $qv (@$qv_ref) {
			if ($qv <= $max) {
				$n++;
			}
		}
	}
	return $n;
}

=head2 contiguous_read_length()

  Usage    : my ($crl_start, $crl_stop) = $ab->contiguous_read_length();
             my ($crl_start, $crl_stop) =
                $ab->contiguous_read_length($windowWidth, $quality_threshold)
  Returns  : the beginning and ending position of the CRL (Contiguous Read Length)

The CRL is (the length of) the longest uninterrupted stretch in a read such that
the average quality of any interval of $windowWidth bases (20 by default) that
is inside such stretch never goes below $threshold (20 by default). The threshold
must be at least 10. The ends of the CRL are further trimmed until there are no bases
with quality values less than 10 within the first five and the last five bases.
The positions are counted from zero. If there is more than one CRL, the position
of the first one is reported.

=cut

sub contiguous_read_length {
	my $self = shift;
	my $window_width = 20;
	$window_width = shift if (@_);
	my $qv_threshold = 20;
	$qv_threshold = shift if (@_);
	
	my $qv_ref = $self->quality_values_ref();
	my $N = scalar(@$qv_ref);
	return (0, 0) if ($N < $window_width);
	my $crl = 0; # 1 if we are inside a crl, 0 otherwise
	my $start = 0;
	my $new_start = 0;
	my $stop = 0;
	my $threshold = $window_width * $qv_threshold;
	my $q = 0;
	my $i;
	for ($i = 0; $i < $window_width; $i++) {
		$q += $$qv_ref[$i];
	}
	$crl = 1 if ($q >= $threshold);
	while ($i < $N) {
		$q -= $$qv_ref[$i - $window_width];
		$q += $$qv_ref[$i];
		if ($crl and $q < $threshold) {
			$crl = 0;
			if ($stop - $start < $i - $new_start - 1) {
				$start = $new_start;
				$stop = $i - 1;
			}
		}
		elsif ( (not $crl) and $q >= $threshold) {
			$crl = 1;
			$new_start = $i - $window_width + 1;
		}
		$i++;
	}
	if ($crl and $stop - $start < $N - $new_start - 1) {
		$start = $new_start;
		$stop = $N - 1;
	}
	my $j = 0;
	while ($start + 4 <= $stop and ($j < 5)) { # Trim the beginning
		 if ($$qv_ref[$start + $j] < 10) {
		 	$start += $j + 1;
		 	$j = 0;
		 }
		 else {
		 	$j++;
		 }
	}
	$j = 0;
	while ($start + 4 <= $stop and ($j < 5)) { # Trim the end
		if ($$qv_ref[$stop - $j] < 10) {
			$stop -= ($j + 1);
			$j = 0;
		}
		else {
			$j++;
		}
	}
	if ($stop - $start < 4) {
		for (my $k = $start; $k <= $stop; $k++) {
			if ($$qv_ref[$k] < 10) {
				return (0, 0);
			}
		}
	}
	return ($start, $stop);
}


=head2 length_of_read()

  Usage    : my $LOR = $ab->length_of_read($windowWidth, $quality_threshold);
             my $LOR = $ab->length_of_read($windowWidth, $quality_threshold, $method);
  Returns  : the Length Of Read (LOR) value, computed using a window of
             $windowWidth bases, a threshold equal to $quality_threshold and,
             optionally, according to the specified $method
             (which, currently, is either the string 'SequencingAnalysis' or the string
             'GoodQualityWindows'). If there are less than $windowWidth quality values,
             the returned value is 0;
             return 0 if there are less than $windowWidth bases;
             returns 0 if the information needed to compute such value is missing.

The Length Of Read (LOR) score gives an approximate measure of the
usable range of high-quality or high-accuracy bases determined by quality
values. Such range can be determined in several ways. Two possible procedures
are currently implemented and described below.

The 'SequencingAnalysis' method (used by default) computes the LOR as the widest
range starting and ending with $windowWidth bases whose average quality 
is greater than or equal to $quality_threshold.

The 'GoodQualityWindows' method computes the LOR as the number of intervals
of $windowWidth bases whose average quality is greater than or equal to 
$quality_threshold.

=cut

sub length_of_read {
	my $self = shift;
	my $window = shift;
	my $qv_threshold = shift;
	my $method = 'SequencingAnalysis';
	if (@_) { $method = shift; }
	
	my $qv_ref = $self->quality_values_ref();
	return 0 unless defined $qv_ref;
	my $LOR = 0;
	my $first = 0;
	my $last = 0;
	my $sum = 0;
	#my $avg = 0;
	
	my $threshold = $window * $qv_threshold;
	my $N = scalar(@$qv_ref);
	if ($N < $window) {
		#print STDERR "Not enough bases to compute LOR.\n";
		#print STDERR "At least $window bases are needed.\n";
		return 0;
	}
	
	# Dispatch according to the chosen method
	if ($method eq 'SequencingAnalysis') {
		# Compute the LOR score as the Sequencing Analysis program does
		# Compute the average quality value in the first window
		my $i;
		# Determine the first window with average qv >= $qv_threshold		
		my $start = 0;
		$sum = 0;
		for ($i = 0; $i < $window; $i++) {
			$sum += $$qv_ref[$i];
		}		
		if ($sum < $threshold) {
			do {
				$sum -= $$qv_ref[$i - $window];
				$sum += $$qv_ref[$i];
				$i++;
			} while ($sum < $threshold and $i < $N);
			$start = $i - $window;
		}
		
		# Determine the last window with average qv >= $qv_threshold
		my $stop = $N - 1;
		$sum = 0;
		for ($i = $stop; $i > $stop - $window; $i--) {
			$sum += $$qv_ref[$i];
		}
		if ($sum < $threshold) {
			do {
				$sum -= $$qv_ref[$i + $window];
				$sum += $$qv_ref[$i];
				$i--;
			} while ($sum < $threshold and $i >= 0);
			$stop = $i + $window;
		}
		
		$LOR = ($stop > $start) ? ($stop - $start + 1) : 0;
	} # end if ($method eq 'SequencingAnalysis')
	else { 
		# This method computes the LOR as the number
		# of windows of size $window having average quality value >= $threshold

		# Compute the average quality value in the first window
		$sum = 0;
		my $i; # Points to the right end of the next window to be processed
		for ($i = 0; $i < $window; $i++) {
			$sum += $$qv_ref[$i];
		}
		#$avg = $sum / $window;
		#if ($avg >= $qv_threshold) {
		if ($sum >= $threshold) {
			$LOR++;
		}
		while ($i < $N) {
			# Compute the average of the shifted window
			$sum -= $$qv_ref[$i - $window];
			$sum += $$qv_ref[$i];
			#$avg = $sum / $window;
			#if ($avg >= $qv_threshold) {
			if ($sum >= $threshold) {
				$LOR++;
			}
			$i++;		
		}
	}
	
	return $LOR;
}

########################################################################

=head1 INTERNAL FUNCTIONS

The following methods are meant for internal use only and should not be used
as user functions.

=head2 _ieee_single_prec_float

  Usage    : _ieee_single_prec_float($string_32bits)
  Returns  : the floating number corresponding to the given 32 bit string.

Interprets the 32 bit string in the standard IEEE format:

  <sign (1 bit)><exponent (8 bits)><mantissa (23 bits)>
  
The value is computed as:

  sign * 1.mantissa * 2**(exponent - 127)

=cut

sub _ieee_single_prec_float {
	my $self = shift;
	my $b = shift;
	my $sign = (substr($b, 0, 1) eq '0') ? 1 : -1;
	my $exp  = $self->_bit_string_as_unsigned_integer(substr($b, 1, 8));
	my $m    = $self->_bit_string_as_decimal_fraction('1.' . substr($b, 9, 23));
	return $sign * $m * (2**($exp - 127));
}


=head2 _bit_string_as_unsigned_integer()

  Usage    : my $n = _bit_string_as_unsigned_integer('10000101');
  Returns  : the decimal value of the given bit string.
  
Interprets the given bit string as an unsigned integer,
and returns its decimal representation (e.g., for example
above, returns 133).

=cut

sub _bit_string_as_unsigned_integer {
	my $self = shift;
	my $bit_string = shift;
	my @bit = split('', $bit_string);
	my $f = 0;
	for (my $i = 0; $i < scalar(@bit); $i++) {
		$f += $bit[$i] * 2**($#bit - $i); 
	}
	return $f;	
};

=head2 _bit_string_as_decimal_fraction()

  Usage    : my $f = _bit_string_as_decimal_fraction('010.01000');
  Returns  : the decimal fraction corresponding to the given
             bit string.

Interprets the given bit string as an unsigned fractional value
and returns the corresponding decimal value. For example,
010.0100 is converted into 2.25.

=cut

sub _bit_string_as_decimal_fraction {
	my $self = shift;
	my $bit_string = shift;
	my ($i, $f) = split('\.', $bit_string);
	my @fbit = split('', $f);
	my $value = $self->_bit_string_as_unsigned_integer($i);
	for (my $j = 0; $j < scalar(@fbit); $j++) {
		$value += $fbit[$j] * 2**(-$j-1); 
	}
	return $value;
};



sub _debug {
	my $self = shift;
	confess "usage: thing->_debug(level)" unless @_ == 1;
	my $level = shift;
	if (ref($self))  {
		$self->{"_DEBUG"} = $level; # just myself
	} 
	else {
		$Debugging = $level; # whole class
	}
}

sub DESTROY {
	my $self = shift;
	if ($Debugging || $self->{"_DEBUG"}) {
		carp "Destroying $self " . $self->name;
	}
}

sub END {
	if ($Debugging) {
		print "All ABIF objects are going away now.\n";
	}
}
    
=head1 AUTHOR

Nicola Vitacolonna, C<< <vitacolonna at appliedgenomics.org> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-bio-trace-abif at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Trace-ABIF>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Trace::ABIF

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Trace-ABIF>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Trace-ABIF>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Trace-ABIF>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Trace-ABIF>

=back

=head1 SEE ALSO

See L<http://www.appliedbiosystems.com/support/> for the ABIF format file
specification sheet.

There is an ABI module on CPAN (L<http://search.cpan.org/~malay/ABI-0.01/>).

bioperl-ext also parses ABIF files and other trace formats.

You are welcome at L<http://www.appliedgenomics.org>!

=head1 ACKNOWLEDGEMENTS

Thanks to Simone Scalabrin for many helpful suggestions and
for the first implementation of the length_of_read() method the way Sequencing Analysis
does it!

Some explication about the way Sequencing Analysis computes some
parameters has been found at L<http://keck.med.yale.edu/dnaseq/>.

=head1 COPYRIGHT & LICENSE

Copyright 2006 Nicola Vitacolonna, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

1; # End of Bio::Trace::ABIF
