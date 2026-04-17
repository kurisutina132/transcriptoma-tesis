SHELL := /bin/bash
PROJECT_DIR := /mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma

.PHONY: transport clean-transport

transport:
	cd $(PROJECT_DIR) && Rscript report_transport_pipeline_full.R

clean-transport:
	rm -rf $(PROJECT_DIR)/reports/transport_full
