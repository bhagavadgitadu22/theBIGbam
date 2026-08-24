from thebigbam.plotting.controls.searchable_select import SearchableSelect


def test_filter_search_contract_is_configurable_and_encoded_in_frontend():
    widget = SearchableSelect(server_search=True, min_search_chars=2)

    assert widget.min_search_chars == 2
    assert widget.search_result_nonce == 0
    assert "query.length >= model.min_search_chars" in widget._esm
    assert "Type at least ${model.min_search_chars} characters" in widget._esm
    assert "pendingFilterLoads" in widget._esm
    assert "if (previous) previous([])" in widget._esm
