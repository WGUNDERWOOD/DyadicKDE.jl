all: test format doc rep

.PHONY: test format doc rep clean

test:
	@echo -e "\e[0;35m\033[1mTesting Julia package...\e[0;30m\033[0m"
	@julia --depwarn=yes --project --color=yes -e 'using Pkg; Pkg.test(coverage=true)'

format:
	@echo -e "\e[0;35m\033[1mFormatting Julia package...\e[0;30m\033[0m"
	@julia --project --color=yes -e 'using JuliaFormatter; format(".")'

doc:
	@echo -e "\e[0;35m\033[1mMaking Julia documentation...\e[0;30m\033[0m"
	@cd docs && julia --project=".." --color=yes make.jl

rep:
	@echo -e "\e[0;35m\033[1mRunning replication files...\e[0;30m\033[0m"
	@cd replication && make

clean:
	@echo -e "\e[0;35m\033[1mCleaning up Julia package...\e[0;30m\033[0m"
	@rm -r docs/build
