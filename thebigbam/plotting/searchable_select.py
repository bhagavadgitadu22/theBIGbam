import param
from panel.custom import JSComponent


class SearchableSelect(JSComponent):
    """Searchable select dropdown using Tom Select."""

    value = param.String(default="")
    options = param.List(default=[])
    placeholder = param.String(default="")

    _stylesheets = [
        "https://cdn.jsdelivr.net/npm/tom-select@2.4.1/dist/css/tom-select.css",
        # Override styles to match Bokeh widget look
        """
        .ts-wrapper { width: 100%; }
        .ts-wrapper .ts-control { border: 1px solid #ccc; border-radius: 4px; padding: 4px 8px; min-height: 31px; font-size: 14px; }
        .ts-dropdown { margin-top: 0 !important; }
        """
    ]

    _esm = """
    import TomSelect from "https://esm.sh/tom-select@2.4.1";

    export function render({ model }) {
        const container = document.createElement('div');
        container.style.width = '100%';
        const select = document.createElement('select');
        select.setAttribute('placeholder', model.placeholder);
        container.appendChild(select);

        let allOptions = model.options.map(o => ({value: o, text: o}));

        // Seed the pre-set value so it displays correctly
        const initialOptions = [];
        if (model.value) {
            initialOptions.push({value: model.value, text: model.value});
        }

        const ts = new TomSelect(select, {
            create: false,
            maxOptions: null,
            placeholder: model.placeholder,
            options: allOptions,
            items: model.value ? [model.value] : [],
            onChange: (val) => { model.value = val; }
        });

        model.on('options', () => {
            const currentVal = model.value;
            allOptions = model.options.map(o => ({value: o, text: o}));
            ts.clearOptions();
            ts.addOptions(allOptions);
            // Restore the current value if still valid
            if (currentVal && model.options.includes(currentVal)) {
                ts.setValue(currentVal, true);
            } else {
                ts.setValue('', true);
            }
        });

        model.on('value', () => {
            // Ensure the option exists before setting value
            if (model.value && !ts.options[model.value]) {
                ts.addOption({value: model.value, text: model.value});
            }
            ts.setValue(model.value, true);
        });

        return container;
    }
    """
