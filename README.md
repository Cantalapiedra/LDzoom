# LDzoom
An R script to reduce a hapmap based on LD decay around a specified marker

Usage: /home/carlospc/bin/LDzoom.R [options]

Options:
	--marker=CHARACTER
		Marker to be used a center of the interval and to compute LD around it.

	--hmp=CHARACTER
		Hapmap to compute the LD and zoom into it.

	--chrom=CHARACTER
		Chromosome where --marker is located (both --chrom, --start and --end are used to specify which region to zoom in for those markers with multiple mappings to different loci.

	--start=CHARACTER
		Start position of --marker alignment to the reference (both --chrom, --start and --end are used to specify which region to zoom in for those markers with multiple mappings to different loci.

	--end=CHARACTER
		End position of --marker alignment to the reference (both --chrom, --start and --end are used to specify which region to zoom in for those markers with multiple mappings to different loci.

	--ldthres=CHARACTER
		When the LD from a given position to --marker goes below this --ldthres, the previous position is a end limit of the final interval.

	--interval_max=CHARACTER
		The final interval reported will never exceed --interval_max, even if the LD does not ever go below --ldthres

	--interval_min=CHARACTER
		The final interval span will never be smaller than --interval_min, even if the LD goes below --ldthres closer to the --marker

	--window_size=CHARACTER
		How many adjacent markers will be used to compute each LD value.

	--window_step=CHARACTER
		How many markers will skip one LD window from the previous LD window.

	-h, --help
		Show this help message and exit

Note: marker, hmp, chrom, start and end are mandatory parameters.
marker, chrom, start and end define a marker locus, which must be included in the region defined by markers in hmp.
ldthres value will help define the new region based on LD decay.
