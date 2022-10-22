local M = {}

local ts_utils = require("nvim-treesitter.ts_utils")

local lua_keys = {
	assignment_statement = true,
	function_call = true,
	arguments = true,
	parameters = true,
	condition = true,
}

local js_keys = {
	variable_declarator = true,
	formal_parameters = "statement_block",
	parenthesized_expression = true,

	sibling_block_placement_behaviour = "inside",
}

local python_keys = {
	expression_statement = true,
	parameters = "block",
	parenthesized_expression = "block",
	comparison_operator = "block",

	sibling_block_placement_behaviour = "before",
	placement_default_behaviour = "behind",
}

M.supported_languages = {
	lua = lua_keys,
	typescript = js_keys,
	typescriptreact = js_keys,
	javascript = js_keys,
	javascriptreact = js_keys,
	python = python_keys,
}

local skip_types = {
	comment = true,
	["=>"] = true,
	[":"] = true,
}

local function getSibling(node)
	local type_text_node = node:type()

	local sibling = node:next_sibling()

	if sibling == nil then
		return node
	end

	while true do
		if skip_types[sibling:type()] then
			sibling = sibling:next_sibling()
		else
			break
		end
	end

	local type_text_sibling = sibling:type()

	if M.language_keys[type_text_node] == type_text_sibling then
		return sibling, M.language_keys["sibling_block_placement_behaviour"]
	else
		return node, M.language_keys["placement_default_behaviour"]
	end
end

local function getParent(node, cnt)
	if cnt > 10 or node == nil then
		return nil
	end

	local type_text = node:type()

	if M.language_keys[type_text] then
		local outnode, direction = getSibling(node)
		return outnode, direction
	end

	return getParent(node:parent(), cnt + 1)
end

local function set_print_statement(node, print_text, placement, format) -- hello lll
	local start_row, start_idx, end_row, end_idx = node:range()

	local row_s
	local row_e
	local updated_text = print_text

	if placement == "inside" then
		row_s = start_row
		row_e = start_row + 1

		local current_text = vim.api.nvim_buf_get_lines(0, row_s, row_e, true)[1]
		updated_text = current_text:sub(1, start_idx + 1) .. print_text .. ";" .. current_text:sub(start_idx + 2)
	elseif placement == "before" then
		row_s = start_row
		row_e = start_row + 1

		local current_text = vim.api.nvim_buf_get_lines(0, row_s, row_e, true)[1]
		updated_text = current_text:sub(1, start_idx + 0) .. print_text .. ";" .. current_text:sub(start_idx + 1)
	elseif placement == "behind" then
		row_s = start_row
		row_e = start_row + 1

		local current_text = vim.api.nvim_buf_get_lines(0, row_s, row_e, true)[1]
		updated_text = current_text:sub(1, end_idx + 0) .. ";" .. print_text .. current_text:sub(end_idx + 1)
	else
		row_s = end_row + 1
		row_e = end_row + 1
	end

	vim.api.nvim_buf_set_lines(0, row_s, row_e, true, { updated_text })

	if format then
		print(row_s, row_e)
		vim.lsp.buf.format({ range = {
			["start"] = { row_s + 0, 0 },
			["end"] = { row_e + 2, -1 },
		} })
	end
end

M.insert_print = function()
	local ft = vim.bo.filetype
	M.language_keys = M.supported_languages[ft]
	if M.language_keys == nil then
		print("No language keys found for filetype: " .. ft)
	end

	local node_at_cursor = ts_utils.get_node_at_cursor(0)
	local outp, placement = getParent(node_at_cursor, 0) --[[ or node_at_cursor ]]

	outp = outp or node_at_cursor
	placement = placement or M.language_keys["placement_default_behaviour"]
	set_print_statement(outp, "print('HELLO')", placement, true)
end

--[[ M.insert_print() ]]

--[[ vim.keymap.set("n", "<leader>lg", function() ]]
--[[   R("zippy").insert_print() ]]
--[[ end) ]]

--[[ local r = ts_locals.get_scope_tree(node_at_cursor) ]]
--[[ local closest_scope = r[1] ]]
--[[ print(closest_scope, "closest_scope") ]]
--[[ print(node_at_cursor, "node_at_cursor") ]]

--[[ for k, v in ipairs(r) do ]]
--[[ 	print(k, v) ]]
--[[ 	P(v) ]]
--[[ 	local text_current = vim.treesitter.query.get_node_text(v, 0) ]]
--[[ 	print(text_current) ]]
--[[ 	ts_utils.highlight_node(v, 0, usage_namespace, "ZIPPY") ]]
--[[ 	break ]]
--[[ end ]]

--[[ ts_utils.highlight_node(current_scope, 0, usage_namespace, "ZIPPY") ]]
--[[ P(current_scope) ]]

--[[ local text_current = vim.treesitter.query.get_node_text(current_scope, 0) ]]
--[[ P(text_current) ]]

--[[ local res = ts_locals.get_scopes() ]]
--[[ P(res) ]]

--[[ local ts_statusline = require("nvim-treesitter.statusline") ]]
--[[  local f = require'nvim-treesitter'.statusline({}) ]]
--[[ print(f, "F") ]]

--[[ local usage_namespace = vim.api.nvim_create_namespace("nvim-treesitter-usages") ]]
--[[ local cursor_pos = vim.api.nvim_win_get_cursor(0) ]]
--[[ vim.api.nvim_set_hl(0, "ZIPPY", { default = true, bg = "#000000", fg = "#ffffff" }) ]]

--[[ local text = vim.treesitter.query.get_node_text(node_at_cursor, 0) ]]
--[[ P(text) ]]

--[[ local parent = node_at_cursor:parent() ]]
--[[ local textParent = vim.treesitter.query.get_node_text(parent, 0) ]]
--[[ P(textParent) ]]

--[[ local parent = ts_utils.get_root_for_node(node_at_cursor) ]]
--[[ local textParent = vim.treesitter.query.get_node_text(parent, 0) ]]
--[[ P(textParent) ]]

--[[ ts_utils.highlight_node(node_at_cursor, 0, usage_namespace, "ZIPPY") ]]
--[[ print(node_at_cursor) ]]
--[[ print(ts_utils.get_node_range(node_at_cursor), "node_length") ]]
--[[ local range_of_node = ts_utils.get_node_range(node_at_cursor) ]]
--[[ P(range_of_node) ]]
--[[ ts_utils.goto_node(node_at_cursor, 1, 1) ]]
--[[ case variable ]]
--[[ start at itentifier and iterate until the next assignemt_statement, local decleration ]]

return M
