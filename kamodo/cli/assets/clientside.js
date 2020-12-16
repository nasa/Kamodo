if (!window.dash_clientside) {
     window.dash_clientside = {}
 }

window.dash_clientside.clientside = {

    figure: function (fig_dict, title) {

        if (!fig_dict) {
            throw "Figure data not loaded, aborting update."
        }

        // Copy the fig_data so we can modify it
        // Is this required? Not sure if fig_data is passed by reference or value
        fig_dict_copy = {...fig_dict};

        fig_dict_copy["layout"]["title"] = title;

        return fig_dict_copy

    },

    subplot: function (fig_dict, input_n) {

        if (!fig_dict) {
            throw "Figure data not loaded, aborting update."
        }

        // Copy the fig_data so we can modify it
        // Is this required? Not sure if fig_data is passed by reference or value
        fig_dict_copy = {...fig_dict};

        var nkeys = Object.keys(fig_dict_copy).length
        
        console.log('nkeys/input_n', Math.ceil(nkeys/input_n));
        console.log('nkeys', nkeys);
        rows = input_n; // nkeys = rows*columns
        columns = Math.ceil(nkeys/input_n);
        console.log('rows', rows);
        console.log('columns', columns);
        var layout = {
          grid: {rows: rows, columns: columns, pattern: 'independent'},
        };


        var traces = [];
        
        console.log('figure keys:', Object.keys(fig_dict_copy));
        var fig_count = 1;

        var xaxis = 'x';
        var yaxis = 'y';

        for (var key in fig_dict_copy){
            
            if (fig_count > 1){
                xaxis = 'x'.concat(fig_count);
                yaxis = 'y'.concat(fig_count);
            }

            var t_ = fig_dict_copy[key];
            t_['xaxis'] = xaxis;
            t_['yaxis'] = yaxis;

            traces.push(t_);
            fig_count ++;
        }
        
        fig_out = {data: traces, layout: layout};

        return fig_out;

    },

    subplots: function (fig_dict, var_names) {

        if (!fig_dict) {
            throw "Figure data not loaded, aborting update."
        }

        // Copy the fig_data so we can modify it
        // Is this required? Not sure if fig_data is passed by reference or value
        fig_dict_copy = {...fig_dict};

        var nkeys = var_names.length;

        rows = nkeys; // nkeys = rows*columns
        columns = 1;
        console.log('rows', rows);
        console.log('columns', columns);
        var layout = {
          grid: {rows: rows, columns: columns, pattern: 'independent'},
        };


        var traces = [];
        
        console.log('figure keys:', Object.keys(fig_dict_copy));
        console.log('variables:', var_names)
        var fig_count = 1;

        var xaxis = 'x';
        var yaxis = 'y';

        for (var k_ in var_names){
            console.log(k_, 'in fig_dict:', k_ in fig_dict_copy);
        }

        for (var key in var_names){
            
            
            if (fig_count > 1){
                xaxis = 'x'.concat(fig_count);
                yaxis = 'y'.concat(fig_count);
            }

            var t_ = fig_dict_copy[var_names[key]];
            t_['xaxis'] = xaxis;
            t_['yaxis'] = yaxis;

            traces.push(t_);
            fig_count ++;
        }
        
        fig_out = {data: traces, layout: layout};

        return fig_out;

    },


}