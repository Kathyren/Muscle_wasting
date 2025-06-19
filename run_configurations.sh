#!/bin/bash

CSV_FILE="UseCases/UseCase0/configurations.csv"
SCRIPT="network/main_pipeline.py" 
CONFIG_DIR="UseCases/UseCase0/configurations_YML"
mkdir -p $CONFIG_DIR

# Read header
read -r HEADER_LINE < "$CSV_FILE"
IFS=',' read -r -a HEADERS <<< "$HEADER_LINE"

# Row counter
i=0

tail -n +2 "$CSV_FILE" | while IFS=',' read -r -a VALUES; do
    ((i++))
    CONFIG_PATH="$CONFIG_DIR/config_$i.yml"
    CMD="python $SCRIPT"

    YAML_CONTENT=""

    declare -A coeffs

    for j in "${!HEADERS[@]}"; do
        key="${HEADERS[j]}"
        value="${VALUES[j]}"
        value="$(echo "$value" | xargs | sed 's/^"//;s/"$//')"

        if [ -n "$value" ]; then
            case "$key" in
                network) YAML_CONTENT+="oNetwork: \"$value\"\n" ;;
                DE_data) YAML_CONTENT+="path_DDS_data: \"$value\"\n" ;;
                Tissue_expression_data) YAML_CONTENT+="path_tissue_data: \"$value\"\n" ;;
                Cell_type_data) YAML_CONTENT+="path_cell_type_data: \"$value\"\n" ;;
                Pathway_enrichment) YAML_CONTENT+="path_pathway_file: \"$value\"\n" ;;
                TF_enrichment) YAML_CONTENT+="path_tf_file: \"$value\"\n" ;;
                Cutoff)
                    CUT=($value)
                    YAML_CONTENT+="page_rank_cutoff: [${CUT[*]}]\n"
                    ;;
                coeff_dds|coeff_tissue|coeff_cellular|coeff_pathway_svd|coeff_tf)
                    # Extract base key
                    base_key="${key#coeff_}"
                    coeffs["$base_key"]="$value"
                    ;;
            esac
        fi
    done

    # Add coefficient block
    COEFF_YAML="coefficients: {"
    for key in dds tissue cellular pathway_svd tf; do
        val="${coeffs[$key]:-1}"  # default to 1
        COEFF_YAML+="'$key': $val, "
    done
    COEFF_YAML="${COEFF_YAML%, }"  # remove trailing comma
    COEFF_YAML+="}\n"

    YAML_CONTENT+="$COEFF_YAML"

    # Save YAML
    echo -e "$YAML_CONTENT" > "$CONFIG_PATH"
    echo "[INFO] Saved config: $CONFIG_PATH"

    # Add --config to command
    CMD+=" --config \"$CONFIG_PATH\""

    # CLI arguments (network_name, Cutoff, output, open_cytoscape)
    for j in "${!HEADERS[@]}"; do
        key="${HEADERS[j]}"
        value="${VALUES[j]}"
        value="$(echo "$value" | xargs | sed 's/^"//;s/"$//')"

        if [ -n "$value" ]; then
            if [[ "$key" == "Cutoff" || "$key" == "network_name" || "$key" == "output" ]]; then
                if [ "$key" == "Cutoff" ]; then
                    CMD+=" --$key $value"
                elif [ "$key" == "open_cytoscape" ]; then
                    [[ "$value" =~ ^(true|True|1|yes)$ ]] && CMD+=" --$key True" || CMD+=" --$key False"
                else
                    CMD+=" --$key \"$value\""
                fi
            fi
        fi
    done

    echo "[INFO] Running: $CMD"
    eval $CMD
done

    echo "Running: $CMD"
    eval $CMD
done

