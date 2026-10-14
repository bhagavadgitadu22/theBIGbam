// Execute the production component with minimal DOM and Tom Select adapters.
const fs = require('fs');
const assert = require('assert').strict;
const {source, name, pasted, current} = JSON.parse(fs.readFileSync(0, 'utf8'));
const handlers = new Map();
const model = {
    options: [], value: '', search_query: '', search_result_nonce: 0,
    server_search: true, allow_custom: false, min_search_chars: 0,
    on(name, callback) { handlers.set(name, callback); },
};
const element = () => ({style: {}, setAttribute() {}, appendChild() {}, addEventListener() {}});
const document = {getElementById: () => ({}), createElement: element};
let ts;
class TomSelect {
    constructor(select, config) {
        ts = this;
        this.settings = config;
        this.options = {};
        this.control_input = {value: '', matches: () => false};
        this.results = [];
    }
    // Tom Select 2.4.1 normalizes input this way before load/onType.
    inputValue() { return this.control_input.value.trim(); }
    addOption(option) { this.options[option.value] = option; }
    removeOption(value) { delete this.options[value]; }
    refreshOptions() {}
    setValue(value) { this.value = value; }
}
const component = source
    .replace(/import TomSelect from [^;]+;/, '')
    .replace('export function render', 'function render');
new Function('TomSelect', 'document', component + '\nreturn render;')(TomSelect, document)({model});

ts.control_input.value = pasted;
const query = ts.inputValue();
ts.settings.onType(query);
let callbackResults;
ts.settings.load(query, results => { callbackResults = results; });
const request = JSON.parse(model.search_request);
assert.equal(request.query, query);
assert.equal(query, name);

// Mimic input changing while the server query is in flight, if requested.
if (current !== null) ts.control_input.value = current;
model.options = [name];
handlers.get('options')();
model.search_result_query = request.query;
model.search_result_nonce = request.nonce;
handlers.get('search_result_nonce')();

const expected = ts.inputValue() === query ? [{value: name, text: name}] : [];
assert.deepEqual(callbackResults, expected);
if (expected.length) {
    assert.equal(ts.options[name].value, name);
    assert.ok(!ts.settings.render.no_results({input: query}, text => text).includes('Loading'));
    // A repeat search should resolve from the cached accepted response.
    callbackResults = undefined;
    ts.settings.load(query, results => { callbackResults = results; });
    assert.deepEqual(callbackResults, expected);
}
