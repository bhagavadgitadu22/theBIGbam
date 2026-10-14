import json
import shutil
import subprocess
from pathlib import Path

import pytest

from thebigbam.plotting.controls.searchable_select import SearchableSelect, decode_search_request


def test_filter_search_contract_is_configurable_and_encoded_in_frontend():
    widget = SearchableSelect(server_search=True, min_search_chars=2)

    assert widget.min_search_chars == 2
    assert widget.search_result_nonce == 0
    assert widget.scope_nonce == 0
    assert "query.length >= model.min_search_chars" in widget._esm
    assert "Type at least ${model.min_search_chars} characters" in widget._esm
    assert "pendingFilterLoads" in widget._esm
    assert "queryResults.get(query)" in widget._esm
    assert "queryResults.set(query, newOptions)" in widget._esm
    assert "pendingInputQuery = null" in widget._esm
    assert "replaceOptions([], false)" in widget._esm
    assert "thebigbam-autocomplete-spinner" in widget._esm
    assert "queryResults.clear()" in widget._esm
    assert "if (previous) previous.callback([])" in widget._esm
    assert "model.on('scope_nonce'" in widget._esm
    assert "dropdownParent: 'body'" in widget._esm
    assert "body > .ts-dropdown { z-index: 10000; }" in widget._esm
    assert "model.search_request = JSON.stringify" in widget._esm
    assert "model.search_nonce" not in widget._esm
    assert "request.nonce !== model.search_result_nonce" in widget._esm
    assert "create: model.allow_custom" in widget._esm
    assert "if (model.allow_custom) model.value = query" in widget._esm


def test_free_text_autocomplete_mode_is_explicit():
    widget = SearchableSelect(server_search=True, allow_custom=True)

    assert widget.allow_custom
    assert widget.options == []


def test_atomic_search_request_keeps_nonce_and_query_together():
    assert decode_search_request('{"nonce": 7, "query": "337R"}') == (7, "337R")
    with pytest.raises(ValueError, match="invalid autocomplete request"):
        decode_search_request("not-json")


@pytest.mark.parametrize("padding", ["", " ", "\t", "\n", "\u00a0"])
@pytest.mark.parametrize("current", [None, "different_contig", ""])
def test_pasted_exact_name_response_uses_normalized_active_query(padding, current):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to execute the autocomplete response handler")
    name = "Ferret_FED_2019_07_19_NODE_485927_length_336492_cov_115.9959_ID_971853"
    result = subprocess.run(
        [node, str(Path(__file__).with_name("searchable_select_harness.cjs"))],
        input=json.dumps({
            "source": SearchableSelect._esm,
            "name": name,
            "pasted": f"{padding}{name}{padding}",
            "current": current,
        }),
        text=True,
        capture_output=True,
        timeout=10,
    )
    assert result.returncode == 0, result.stderr
