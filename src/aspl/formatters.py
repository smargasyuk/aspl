from .events import CassetteExon, SpliceJunction, SpliceSite, Exon
from .events import Strand, CoordinateSystem, SpliceSiteType
from .transcripts import Transcript
from dataclasses import dataclass
from typing import TypeVar, cast
import itertools

T = TypeVar("T", CassetteExon, SpliceJunction, SpliceSite, Exon)


@dataclass
class UnderscoreSeparated:
    coordinate_system: CoordinateSystem = CoordinateSystem.ONE_BASED
    include_site_type: bool = False


    def format(self, value: CassetteExon | SpliceJunction | SpliceSite | Exon | Transcript) -> str:
        assert value.is_valid()
        sites = value.get_sites_sorted_by_coordinate()
        coords = [site.coord for site in sites]
        # add proper checks?
        is_splice_site = len(coords) == 1

        if is_splice_site and self.coordinate_system == CoordinateSystem.ZERO_BASED:
            coords[0] = coords[0] - 1
        fields = [sites[0].seqname] + [str(c) for c in coords] + [sites[0].strand]
        if is_splice_site and self.include_site_type:
            fields += [sites[0].type]
        return "_".join(fields)


# seqname should not contain any underscores, needed for transcripts
    def parse(self, input: str, target_type: type[T]) -> T:
        match target_type:

            case t if t is SpliceSite:
                if self.include_site_type:
                    seqname, coord, strand, site_type = input.rsplit("_", maxsplit=3)
                    site_type = SpliceSiteType(site_type)
                else:
                    seqname, coord, strand = input.rsplit("_", maxsplit=2)
                    site_type = SpliceSiteType.UNKNOWN

                coord = int(coord)
                if self.coordinate_system == CoordinateSystem.ZERO_BASED:
                    coord += 1

                strand = Strand(strand)
                output = SpliceSite(seqname, coord, strand, site_type)
                return cast(T, output)


            case t:
                fields = input.split('_')
                strand = Strand(fields[-1])
                seqname = fields[0]
                coords = [int(f) for f in fields[1:-1]]
                if strand == Strand.MINUS:
                    coords = coords[::-1]
                if t is Transcript:
                    splice_site_types = itertools.cycle([SpliceSiteType.ACCEPTOR, SpliceSiteType.DONOR])
                else:
                    splice_site_types = t.splice_site_types
                sites = [SpliceSite(seqname, coord, strand, sst) for coord, sst in zip(coords, splice_site_types)]
                output =  t.from_splice_sites(sites)
                return cast(T, output)  


            # case _:
            #     raise TypeError(f"Unsupported type: {target_type}")


@dataclass
class Bed12Formatter:
    intron_flanks_lenght: int = 3
    BED12TEMPLATE = "{chrom}\t{chromStart}\t{chromEnd}\t{name}\t{score}\t{strand}\t{thickStart}\t{thickEnd}\t{itemRgb}\t{blockCount}\t{blockSizes}\t{blockStarts}"
    name_formatter = UnderscoreSeparated()


    def format(self, value: CassetteExon | SpliceJunction | SpliceSite | Exon, name: str = None, color: str = "0,0,0") -> str:
        match value:
            case CassetteExon():
                name = name if name is not None else self.name_formatter.format(value)
                sites = value.get_sites_sorted_by_coordinate()
                start_coord = sites[0].coord - self.intron_flanks_lenght
                end_coord = sites[-1].coord + self.intron_flanks_lenght - 1
                return self.BED12TEMPLATE.format(
                    chrom = sites[0].seqname,
                    chromStart=start_coord,
                    chromEnd=end_coord,
                    name=name,
                    score="0",
                    strand=sites[0].strand,
                    thickStart=start_coord,
                    thickEnd=end_coord,
                    itemRgb=color,
                    blockCount=3,
                    blockSizes=f"{self.intron_flanks_lenght},{sites[2].coord - sites[1].coord + 1},{self.intron_flanks_lenght}",
                    blockStarts=f"0,{sites[1].coord - sites[0].coord + self.intron_flanks_lenght - 1},{sites[3].coord - sites[0].coord + self.intron_flanks_lenght - 1}"
                )
            case SpliceJunction():
                name = name if name is not None else self.name_formatter.format(value)
                sites = value.get_sites_sorted_by_coordinate()
                start_coord = sites[0].coord - self.intron_flanks_lenght
                end_coord = sites[-1].coord + self.intron_flanks_lenght - 1
                return self.BED12TEMPLATE.format(
                    chrom = sites[0].seqname,
                    chromStart=start_coord,
                    chromEnd=end_coord,
                    name=name,
                    score="0",
                    strand=sites[0].strand,
                    thickStart=start_coord,
                    thickEnd=end_coord,
                    itemRgb=color,
                    blockCount=2,
                    blockSizes=f"{self.intron_flanks_lenght},{self.intron_flanks_lenght}",
                    blockStarts=f"0,{sites[-1].coord - sites[0].coord + self.intron_flanks_lenght - 1}"
                )     
            case SpliceSite():
                name = name if name is not None else self.name_formatter.format(value)
                start_coord = value.coord - 1
                end_coord = value.coord
                return self.BED12TEMPLATE.format(
                    chrom = value.seqname,
                    chromStart=start_coord,
                    chromEnd=end_coord,
                    name=name,
                    score="0",
                    strand=value.strand,
                    thickStart=start_coord,
                    thickEnd=end_coord,
                    itemRgb=color,
                    blockCount=1,
                    blockSizes=f"1",
                    blockStarts=f"0"
                )
            case Exon():
                name = name if name is not None else self.name_formatter.format(value)
                sites = value.get_sites_sorted_by_coordinate()
                start_coord = sites[0].coord - 1
                end_coord = sites[-1].coord
                return self.BED12TEMPLATE.format(
                    chrom = sites[0].seqname,
                    chromStart=start_coord,
                    chromEnd=end_coord,
                    name=name,
                    score="0",
                    strand=sites[0].strand,
                    thickStart=start_coord,
                    thickEnd=end_coord,
                    itemRgb=color,
                    blockCount=1,
                    blockSizes=f"{end_coord-start_coord}",
                    blockStarts=f"0"
                )                                           
            case _:
                raise NotImplementedError
