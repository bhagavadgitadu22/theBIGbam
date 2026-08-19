from thebigbam.plotting.application.availability_titles import AvailabilityCounts, title_html


def test_titles_are_formatted_from_plain_counts():
    titles = title_html(AvailabilityCounts(12, 3, 25, 4, 7), has_mags=True)
    assert "12 available" in titles["contig"]
    assert "3 available" in titles["sample"]
    assert "25 contig/sample pairs, 7 MAG/sample pairs" in titles["filtering"]
    assert "4 available" in titles["mag"]


def test_non_mag_titles_omit_mag_counts():
    titles = title_html(AvailabilityCounts(2, 1, 4), has_mags=False)
    assert "mag" not in titles
    assert "MAG/sample" not in titles["filtering"]
